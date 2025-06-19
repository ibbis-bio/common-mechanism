#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Script that checks the output from hmmscan and prints to screen the results

Usage:
    python check_benign.py -i INPUT -s SEQUENCE -d DATABASE FOLDER
      -i, --input = input sample name (will check for sample.benign.hmmscan file)
      -s, --sequence = input sequence file
      -d, --database = database folder location/path (will check for benign_annotations.tsv)
"""
import logging
import pandas as pd
from commec.tools.blast_tools import get_top_hits
from commec.tools.blastn import BlastNHandler  # For has_hits.
from commec.tools.hmmer import HmmerHandler
from commec.tools.cmscan import CmscanHandler
from commec.config.query import Query
from commec.tools.hmmer import (
    recalculate_hmmer_query_coordinates,
    append_nt_querylength_info
)
from commec.config.result import (
    HitResult,
    ScreenStep,
    ScreenStatus,
    HitScreenStatus,
    MatchRange,
    compare
)

# Constants determining Commec's sensitivity for benign screen.
BENIGN_PROTEIN_EVALUE_CUTOFF : float = 1e-20
MINIMUM_PEPTIDE_COVERAGE : int = 50 # Number is counted in NTs, not AA's.
MINIMUM_QUERY_COVERAGE_FRACTION : float = 0.80
MINIMUM_RNA_BASEPAIR_COVERAGE : int = 50
MINIMUM_SYNBIO_COVERAGE_FRACTION : float = 0.80

logger = logging.getLogger(__name__)

def _filter_benign_proteins(query : Query,
                            hit : HitResult,
                            region : MatchRange,
                            benign_protein_for_query : pd.DataFrame,
                            benign_descriptions : pd.DataFrame) -> list:
    """
    Performs the logic of a benign housekeeping protein dataframe for a
    specific query, and determines whether it contains the criteria to clear 
    a specified hit, for a specified region.

    Benign proteins should have a coverage of at least 80% of the region they are 
    aiming to cover, as well as be at least 50 nucleotide in length.
    """
    logger.debug("\t\tChecking hit %s for benign proteins: ", hit.name)

    # Ignore this region, if there are no overlapping hits.
    benign_protein_for_query_trimmed = _trim_to_region(benign_protein_for_query, region).copy()
    if benign_protein_for_query_trimmed.empty:
        logger.debug("\t\tNo overlapping benign Protein regions for %s", hit.name)
        return []

    benign_protein_for_query_trimmed = _calculate_coverage(benign_protein_for_query_trimmed, region)
    logger.debug("\tCoverage Data for %s: shape: %s preview:\n%s", 
                     hit.name, benign_protein_for_query_trimmed.shape, 
                     benign_protein_for_query_trimmed[["coverage_nt", "coverage_ratio"]].head())
    # Ensure that this benign hit covers at least 50 nucleotides of the query
    benign_protein_for_query_trimmed = benign_protein_for_query_trimmed[
        benign_protein_for_query_trimmed["coverage_nt"] > MINIMUM_PEPTIDE_COVERAGE]

    if benign_protein_for_query_trimmed.empty:
        logger.info("\t --> Insufficient coverage (<%i bp) for housekeeping protein to clear %s (%i-%i).",
        int(MINIMUM_PEPTIDE_COVERAGE), hit.name, region.query_start, region.query_end)
        return []

    # Ensure that this benign hit covers at least 80% of the regulated region.
    benign_protein_for_query_trimmed = benign_protein_for_query_trimmed[
        benign_protein_for_query_trimmed["coverage_ratio"] > MINIMUM_QUERY_COVERAGE_FRACTION]

    benign_protein_for_query_trimmed = benign_protein_for_query_trimmed.reset_index(drop=True)

    if benign_protein_for_query_trimmed.empty:
        logger.info("\t --> Insufficient coverage (<%i%%) for housekeeping protein to clear %s (%i-%i).",
        int(MINIMUM_QUERY_COVERAGE_FRACTION*100), hit.name, region.query_start, region.query_end)
        return []

    # Report top hit for Protein / RNA / Synbio
    benign_hit = benign_protein_for_query_trimmed["subject title"][0]
    benign_hit_description = str(*benign_descriptions["Description"][benign_descriptions["ID"] == benign_hit])
    match_ranges = [
        MatchRange(
        float(benign_protein_for_query_trimmed['evalue'][0]),
        int(benign_protein_for_query_trimmed['s. start'][0]), int(benign_protein_for_query_trimmed['s. end'][0]),
        int(benign_protein_for_query_trimmed['q. start'][0]), int(benign_protein_for_query_trimmed['q. end'][0])
        )
    ]
    benign_hit_outcome = HitResult(
            HitScreenStatus(
                ScreenStatus.PASS,
                ScreenStep.BENIGN_PROTEIN
            ),
            benign_hit,
            benign_hit_description,
            match_ranges,
            annotations={"Coverage: ":float(benign_protein_for_query_trimmed['coverage_ratio'][0])}
        )

    # Rarely, something can be cleared that is already cleared, no need to report on that.
    if hit.recommendation.status not in {ScreenStatus.CLEARED_FLAG, ScreenStatus.CLEARED_WARN}:
        logger.info("\t --> Clearing %s (%s) with house-keeping protein, %s",
                        hit.name, hit.recommendation.status, benign_hit_outcome.name)
        hit.recommendation.status = hit.recommendation.status.clear()
        benign_hit_outcome.recommendation.status = hit.recommendation.status

    return [benign_hit_outcome]

def _filter_benign_rna(query : Query,
                            hit : HitResult,
                            region : MatchRange,
                            benign_rna_for_query : pd.DataFrame) -> list:
    logger.debug("\t\t\tChecking query (%s) hit %s for benign RNA",  query.name, hit.name)
    
    # Filter benign RNA for relevance...
    benign_rna_for_query_trimmed = _trim_to_region(benign_rna_for_query, region).copy()
    if benign_rna_for_query_trimmed.empty:
        logger.debug("\t\t\tNo overlapping benign RNA regions for %s (%i-%i)",
                     hit.name, region.query_start, region.query_end)
        return []

    benign_rna_for_query_trimmed = _calculate_coverage(benign_rna_for_query_trimmed, region)
    logger.debug("\t\t\tCoverage Data for %s: shape: %s preview:\n%s", 
                     hit.name, benign_rna_for_query_trimmed.shape, 
                     benign_rna_for_query_trimmed[["coverage_nt", "coverage_ratio"]].head())
    
    benign_rna_for_query_passed = benign_rna_for_query_trimmed[
        (region.length() - benign_rna_for_query_trimmed["coverage_nt"]) < MINIMUM_RNA_BASEPAIR_COVERAGE]
    benign_rna_for_query_passed = benign_rna_for_query_passed.reset_index(drop=True)

    logger.debug("\t\tSummary Benign RNA for %s: shape: %s preview:\n%s", 
                hit.name, benign_rna_for_query_trimmed.shape,
                benign_rna_for_query_trimmed[["coverage_nt", "coverage_ratio"]].head())
    
    if not benign_rna_for_query_passed.empty:
        benign_hit = benign_rna_for_query_trimmed["subject title"][0]
        benign_hit_description =  benign_rna_for_query_trimmed["description of target"][0]
        match_ranges = [
            MatchRange(
            float(benign_rna_for_query_trimmed['evalue'][0]),
            int(benign_rna_for_query_trimmed['s. start'][0]), int(benign_rna_for_query_trimmed['s. end'][0]),
            int(benign_rna_for_query_trimmed['q. start'][0]), int(benign_rna_for_query_trimmed['q. end'][0])
            )
        ]
        benign_hit_outcome = HitResult(
                HitScreenStatus(
                    ScreenStatus.PASS,
                    ScreenStep.BENIGN_RNA
                ),
                benign_hit,
                benign_hit_description,
                match_ranges,
            )
        
        if hit.recommendation.status not in {ScreenStatus.CLEARED_FLAG, ScreenStatus.CLEARED_WARN}:
            logger.info("\t --> Clearing %s %s (region %i-%i), with Benign RNA %s",
                        hit.recommendation.status, hit.name, region.query_start, region.query_end, benign_hit_outcome.name)
            hit.recommendation.status = hit.recommendation.status.clear()
            benign_hit_outcome.recommendation.status = hit.recommendation.status

        return [benign_hit_outcome]
    
    logger.info("Clear failed for %s (%s) as Benign RNA >%i bases unaccounted for.",
                    hit.name, hit.recommendation.status, MINIMUM_RNA_BASEPAIR_COVERAGE)
    return []

def _filter_benign_dna(query : Query,
                          hit : HitResult,
                          region : MatchRange,
                          benign_dna_for_query : pd.DataFrame) -> list:
    logger.debug("\t\t\tChecking query (%s) hit %s for benign Synbio",  query.name, hit.name)

    # Filter benign SynBio for relevance...
    benign_dna_for_query_trimmed = _trim_to_region(benign_dna_for_query, region).copy()
    if benign_dna_for_query_trimmed.empty:
        logger.debug("\t ... No overlapping benign DNA regions for %s (%i-%i)",
                     hit.name, region.query_start, region.query_end)
        return []

    benign_dna_for_query_trimmed = _calculate_coverage(benign_dna_for_query_trimmed, region)
    logger.debug("\t\t\tCoverage Data for %s: shape: %s preview:\n%s", 
                     hit.name, benign_dna_for_query_trimmed.shape, 
                     benign_dna_for_query_trimmed[["coverage_nt", "coverage_ratio"]].head())

    benign_dna_for_query_trimmed = benign_dna_for_query_trimmed[
        benign_dna_for_query_trimmed["coverage_ratio"] > MINIMUM_SYNBIO_COVERAGE_FRACTION]
    benign_dna_for_query_trimmed = benign_dna_for_query_trimmed.reset_index(drop=True)
    logger.debug("\t\t\tFiltered Coverage Synbio for %s: shape: %s preview:\n%s",
                     hit.name, benign_dna_for_query_trimmed.shape,
                     benign_dna_for_query_trimmed.head())

    if benign_dna_for_query_trimmed.empty:
        logger.info("\t --> Synbio sequences <80%% coverage achieved over hit %s over %i-%i.",
                        hit.name, region.query_start, region.query_end)
        return []

    benign_hit = benign_dna_for_query_trimmed["subject title"][0]
    benign_hit_description =  benign_dna_for_query_trimmed["subject title"][0]
    match_ranges = [
        MatchRange(
        float(benign_dna_for_query_trimmed['evalue'][0]),
        int(benign_dna_for_query_trimmed['s. start'][0]), int(benign_dna_for_query_trimmed['s. end'][0]),
        int(benign_dna_for_query_trimmed['q. start'][0]), int(benign_dna_for_query_trimmed['q. end'][0])
        )
    ]
    benign_hit_outcome = HitResult(
            HitScreenStatus(
                ScreenStatus.PASS,
                ScreenStep.BENIGN_DNA
            ),
            benign_hit,
            benign_hit_description,
            match_ranges,
        )
    
    logger.debug("Processing Benign Hit: %s", benign_hit_outcome)
    
    if hit.recommendation.status not in {ScreenStatus.CLEARED_FLAG, ScreenStatus.CLEARED_WARN}:
        logger.info("\t --> Clearing %s %s region %i-%i as Benign DNA, with synthetic biology part %s",
                        hit.recommendation.status, hit.name, region.query_start, region.query_end, benign_hit_outcome.name)
        hit.recommendation.status = hit.recommendation.status.clear()
        benign_hit_outcome.recommendation.status = hit.recommendation.status
    return [benign_hit_outcome]

def _update_benign_data_for_query(query : Query,
                                  benign_protein : pd.DataFrame,
                                  benign_rna : pd.DataFrame,
                                  benign_dna : pd.DataFrame,
                                  benign_descriptions : pd.DataFrame):
    """
    For a single query, look at all three benign database outputs, and update the 
    single queries hit descriptions to record all the benign hits, as well as clear
    any overlapping WARN or FLAG hits.
    """

    logger.debug("\t\tProcessing benign data for query: %s", query.name)

    benign_protein_for_query = pd.DataFrame()
    benign_rna_for_query = pd.DataFrame()
    benign_dna_for_query = pd.DataFrame()

    logger.info(" Benign checks for %s", query.name)

    # We only care about the benign data for this query, but we need to check for size.
    if not benign_protein.empty:
        benign_protein_for_query = benign_protein[
            benign_protein["query name"].str.rsplit("_", n=1).str[0] == query.name
        ]

    if not benign_rna.empty:
        benign_rna_for_query = benign_rna[
            benign_rna["query name"] == query.name
        ]

    if not benign_dna.empty:
        benign_dna_for_query = benign_dna[
            benign_dna["query acc."] == query.name
        ]

    new_benign_hits = []
    # Separated for logging purposes only.
    new_benign_protein_hits = []
    new_benign_rna_hits = []
    new_benign_dna_hits = []

    # Check every region, of every hit that is a FLAG or WARN, against the Benign screen outcomes.
    for hit in query.result_handle.hits.values():
        # Ignore regions that don't require clearing...
        if ((hit.recommendation.status not in {
            ScreenStatus.FLAG,
            ScreenStatus.WARN
            })
            or
            (hit.recommendation.from_step == ScreenStep.BIORISK)
            ):
            logger.debug("\t\t\tIgnoring %s [%s], not need to clear with benign.", hit.name, hit.recommendation.status)
            continue

        for region in hit.ranges:

            if not benign_protein_for_query.empty:
                query.confirm_has_hits()
                new_benign_protein_hits.extend(
                    _filter_benign_proteins(query, hit, region,
                                            benign_protein_for_query,
                                            benign_descriptions)
                    )
            
            if not benign_rna_for_query.empty:
                query.confirm_has_hits()
                new_benign_rna_hits.extend(
                    _filter_benign_rna(query, hit, region,
                                       benign_rna_for_query)
                )
                
            if not benign_dna_for_query.empty:
                query.confirm_has_hits()
                new_benign_dna_hits.extend(
                    _filter_benign_dna(query, hit, region,
                                          benign_dna_for_query)
                    )
                
    # Report Query Specific lack of hits.
    if len(new_benign_protein_hits) == 0:
        logger.info("\t ... no housekeeping protein hits.")
    if len(new_benign_rna_hits) == 0:
        logger.info("\t ... no benign RNA hits.")
    if len(new_benign_dna_hits) == 0:
        logger.info("\t ... no benign DNA hits.")

    new_benign_hits.extend(new_benign_protein_hits)
    new_benign_hits.extend(new_benign_rna_hits)
    new_benign_hits.extend(new_benign_dna_hits)

    logger.debug("\tNew benign hits added: %i", len(new_benign_hits))
    for benign_addition in new_benign_hits:
        query.result_handle.add_new_hit_information(benign_addition)
        logger.debug("\t\tAdding Benign Hit: %s", benign_addition)

def update_benign_data_from_database(benign_protein_handle : HmmerHandler,
                                     benign_rna_handle : CmscanHandler,
                                     benign_dna_handle : BlastNHandler,
                                     queries : dict[str,Query],
                                     benign_desc : pd.DataFrame):
    """
    Parse the outputs from the protein, rna, and synbio database searches, and populate
    the benign hits into a Screen dataset. Marks those hits that are cleared for benign
    as cleared if benign screen passes them.
    """
    # Reading empty outcomes should result in empty DataFrames, not errors.
    benign_protein_screen_data = benign_protein_handle.read_output()
    append_nt_querylength_info(benign_protein_screen_data, queries)
    recalculate_hmmer_query_coordinates(benign_protein_screen_data)
    logger.debug("\tBenign Protein Data: shape: %s preview:\n%s",
                 benign_protein_screen_data.shape, benign_protein_screen_data.head())
    
    benign_rna_screen_data = benign_rna_handle.read_output()
    logger.debug("\tBenign RNA Data: shape: %s preview:\n%s",
                benign_rna_screen_data.shape, benign_rna_screen_data.head())
    
    benign_dna_screen_data = benign_dna_handle.read_output()
    benign_dna_screen_data = get_top_hits(benign_dna_screen_data)
    logger.debug("\tBenign Synbio Top Hits Data: shape: %s preview:\n%s",
                benign_dna_screen_data.shape, benign_dna_screen_data.head())
    
    for query in queries.values():

        skip = True
        for hit in query.result_handle.hits.values():
            if hit.recommendation.status in {
                ScreenStatus.FLAG,
                ScreenStatus.WARN
            }:
                skip = False

        if skip:
            query.result_handle.recommendation.benign_status = ScreenStatus.SKIP
        
        if query.result_handle.recommendation.benign_status == ScreenStatus.SKIP:
            logger.debug("Skipping query %s, no regulated regions to clear.", query.name)
            continue
        
        _update_benign_data_for_query(query,
                                      benign_protein_screen_data,
                                      benign_rna_screen_data,
                                      benign_dna_screen_data,
                                      benign_desc)

        # Calculate the Benign Screen outcomes for each query.
        query.result_handle.recommendation.benign_status = ScreenStatus.PASS
        # If any hits are still warnings, or flags, propagate that to the benign step.
        for flagged_hit in query.result_handle.get_flagged_hits():
            query.result_handle.recommendation.benign_status = compare(
                flagged_hit.recommendation.status,
                query.result_handle.recommendation.benign_status
                )

def _trim_to_region(data : pd.DataFrame, region : MatchRange):
    datatrim = data[
        ~(
            (data["q. start"] > region.query_end)
            & (data["q. end"] > region.query_end)
        )
        & ~(
            (data["q. start"] < region.query_start)
            & (data["q. end"] < region.query_start)
        )
    ]
    return datatrim

def _calculate_coverage(data : pd.DataFrame, region : MatchRange):
    """ 
    Mutates data to add coverage compared to a specific region, 
    as a ratio and absolute number of nucleotides.
    """
    data.loc[:,"coverage_nt"] = (
        data["q. end"].clip(upper=region.query_end) - 
        data["q. start"].clip(lower=region.query_start)
        )
    data.loc[:,"coverage_ratio"] = data["coverage_nt"] / region.length()
    return data