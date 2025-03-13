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
    # Ignore this region, if there are no overlapping hits.
    benign_protein_for_query_trimmed = _trim_to_region(benign_protein_for_query, region).copy()
    if benign_protein_for_query_trimmed.empty:
        return []

    benign_protein_for_query_trimmed = _calculate_coverage(benign_protein_for_query_trimmed, region)

    # Ensure that this benign hit covers at least 50 nucleotides of the query
    benign_protein_for_query_trimmed = benign_protein_for_query_trimmed[
        benign_protein_for_query_trimmed["coverage_nt"] > MINIMUM_PEPTIDE_COVERAGE]

    if benign_protein_for_query_trimmed.empty:
        logging.info("Not enough housekeeping protein coverage to clear %s found in %s",
        hit.name, query.name)
        return []

    # Ensure that this benign hit covers at least 80% of the regulated region.
    benign_protein_for_query_trimmed = benign_protein_for_query_trimmed[
        benign_protein_for_query_trimmed["coverage_ratio"] > MINIMUM_QUERY_COVERAGE_FRACTION]

    benign_protein_for_query_trimmed = benign_protein_for_query_trimmed.reset_index(drop=True)

    if benign_protein_for_query_trimmed.empty:
        logging.info("Housekeeping protein to clear %s, insufficiently covers query %s",
        hit.name, query.name)
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
        logging.info("Clearing %s (%s) as house-keeping protein, for %s", 
                        hit.name, hit.recommendation.status, query.name)
        hit.recommendation.status = hit.recommendation.status.clear()

    return [benign_hit_outcome]



def _filter_benign_rna(query : Query,
                            hit : HitResult,
                            region : MatchRange,
                            benign_rna_for_query : pd.DataFrame) -> list:
    # Filter benign RNA for relevance...
    benign_rna_for_query_trimmed = _trim_to_region(benign_rna_for_query, region).copy()
    if benign_rna_for_query_trimmed.empty:
        return []

    benign_rna_for_query_trimmed = _calculate_coverage(benign_rna_for_query_trimmed, region)

    #benign_rna_for_query_trimmed = benign_rna_for_query_trimmed.assign(
    #            coverage=region.length() - abs(benign_rna_for_query_trimmed["q. end"] - benign_rna_for_query_trimmed["q. start"])
    #        )
    
    benign_rna_for_query_passed = benign_rna_for_query_trimmed[
        (region.length() - benign_rna_for_query_trimmed["coverage_nt"]) < MINIMUM_RNA_BASEPAIR_COVERAGE]
    benign_rna_for_query_passed = benign_rna_for_query_passed.reset_index(drop=True)

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
            logging.info("Clearing %s (%s) as <%i bases unaccounted for Benign RNA, for %s",
                        hit.name, hit.recommendation.status, MINIMUM_RNA_BASEPAIR_COVERAGE, query.name)
            hit.recommendation.status = hit.recommendation.status.clear()
        return [benign_hit_outcome]
    
    logging.info("Clear failed for %s (%s) as Benign RNA >%i bases unaccounted for, for %s", 
                    hit.name, hit.recommendation.status, MINIMUM_RNA_BASEPAIR_COVERAGE, query.name)
    return []

def _filter_benign_synbio(query : Query,
                          hit : HitResult,
                          region : MatchRange,
                          benign_synbio_for_query : pd.DataFrame) -> list:

    # Filter benign SynBio for relevance...
    benign_synbio_for_query_trimmed = _trim_to_region(benign_synbio_for_query, region).copy()
    if benign_synbio_for_query_trimmed.empty:
        return []
    

    benign_synbio_for_query_trimmed = _calculate_coverage(benign_synbio_for_query_trimmed, region)

    # This original way of doing things is coverage for the entire query.
    #benign_synbio_for_query_trimmed = benign_synbio_for_query_trimmed[
    #    benign_synbio_for_query_trimmed["q. coverage"] > MINIMUM_SYNBIO_COVERAGE_FRACTION]
    #benign_synbio_for_query_trimmed = benign_synbio_for_query_trimmed.reset_index(drop=True)

    benign_synbio_for_query_trimmed = benign_synbio_for_query_trimmed[
        benign_synbio_for_query_trimmed["coverage_ratio"] > MINIMUM_SYNBIO_COVERAGE_FRACTION]
    benign_synbio_for_query_trimmed = benign_synbio_for_query_trimmed.reset_index(drop=True)


    if benign_synbio_for_query_trimmed.empty:
        logging.info("Synbio sequences <80%% coverage achieved over hit %s for query %s.", hit.name, query.name)
        return []

    benign_hit = benign_synbio_for_query_trimmed["subject title"][0]
    benign_hit_description =  benign_synbio_for_query_trimmed["subject title"][0]
    match_ranges = [
        MatchRange(
        float(benign_synbio_for_query_trimmed['evalue'][0]),
        int(benign_synbio_for_query_trimmed['s. start'][0]), int(benign_synbio_for_query_trimmed['s. end'][0]),
        int(benign_synbio_for_query_trimmed['q. start'][0]), int(benign_synbio_for_query_trimmed['q. end'][0])
        )
    ]
    benign_hit_outcome = HitResult(
            HitScreenStatus(
                ScreenStatus.PASS,
                ScreenStep.BENIGN_SYNBIO
            ),
            benign_hit,
            benign_hit_description,
            match_ranges,
        )
    
    if hit.recommendation.status not in {ScreenStatus.CLEARED_FLAG, ScreenStatus.CLEARED_WARN}:
        logging.info("Clearing %s (%s) as Synthetic, for %s",
                        hit.name, hit.recommendation.status, query.name)
        hit.recommendation.status = hit.recommendation.status.clear()

    return [benign_hit_outcome]

def _update_benign_data_for_query(query : Query,
                                  benign_protein : pd.DataFrame,
                                  benign_rna : pd.DataFrame,
                                  benign_synbio : pd.DataFrame,
                                  benign_descriptions : pd.DataFrame):
    """
    For a single query, look at all three benign database outputs, and update the 
    single queries hit descriptions to record all the benign hits, as well as clear
    any overlapping WARN or FLAG hits.
    """
    # We only care about the benign data for this query.
    benign_protein_for_query = benign_protein[
        benign_protein["query name"].str.rsplit("_", n=1).str[0] == query.name
    ]
    benign_rna_for_query = benign_rna[
        benign_rna["query name"] == query.name
    ]
    benign_synbio_for_query = benign_synbio[
        benign_synbio["query acc."] == query.name
    ]

    new_benign_hits = []

    # Early report if query has no benign hits per category:
    if benign_protein_for_query.empty:
        logging.info("\t...no housekeeping protein hits for %s", query.name)
    if benign_rna_for_query.empty:
        logging.info("\t...no benign RNA hits for %s", query.name)
    if benign_synbio_for_query.empty:
        logging.info("\t...no Synbio sequence hits for %s", query.name)

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
            continue
        for region in hit.ranges:

            if not benign_protein_for_query.empty:
                new_benign_hits.extend(
                    _filter_benign_proteins(query, hit, region,
                                            benign_protein_for_query,
                                            benign_descriptions)
                    )
            
            if not benign_rna_for_query.empty:
                new_benign_hits.extend(
                    _filter_benign_rna(query, hit, region,
                                       benign_rna_for_query)
                )
                
            if not benign_synbio_for_query.empty:
                new_benign_hits.extend(
                    _filter_benign_synbio(query, hit, region,
                                          benign_synbio_for_query)
                    )

    # We cannot alter the hits dictionary whilst iterating,
    # So we add everything afterwards.
    for benign_addition in new_benign_hits:
        query.result_handle.add_new_hit_information(benign_addition)


def update_benign_data_from_database(benign_protein_handle : HmmerHandler,
                                     benign_rna_handle : CmscanHandler,
                                     benign_synbio_handle : BlastNHandler,
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

    benign_rna_screen_data = benign_rna_handle.read_output()
    benign_synbio_screen_data = benign_synbio_handle.read_output()
    benign_synbio_screen_data = get_top_hits(benign_synbio_screen_data)

    for query in queries.values():
        _update_benign_data_for_query(query,
                                      benign_protein_screen_data,
                                      benign_rna_screen_data,
                                      benign_synbio_screen_data,
                                      benign_desc)

        # Calculate the Benign Screen outcomes for each query.
        query.result_handle.recommendation.benign_status = ScreenStatus.PASS
        # If any hits are still warnings, or flags, propagate that.
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
