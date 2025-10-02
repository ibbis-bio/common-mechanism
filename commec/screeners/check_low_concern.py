#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Script that checks the output from 3 databases of low concern, deriving protein, 
dna, and rna sources. Protein (hmmscan), rna (cmscan), and dna (blastn) are 
checked via the entry point `parse_low_concern_hits``.
The input dictionary of Query objects is updated in place with the resulting 
relevant low concern information.
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
)

# Constants determining Commec's sensitivity for low_concern screen.
LOW_CONCERN_PROTEIN_EVALUE_CUTOFF : float = 1e-20
MINIMUM_PEPTIDE_COVERAGE : int = 50 # Number is counted in NTs, not AA's.
MINIMUM_QUERY_COVERAGE_FRACTION : float = 0.80
MINIMUM_RNA_BASEPAIR_COVERAGE : int = 50
MINIMUM_SYNBIO_COVERAGE_FRACTION : float = 0.80

logger = logging.getLogger(__name__)

def _filter_low_concern_proteins(query : Query,
                            hit : HitResult,
                            region : MatchRange,
                            low_concern_protein_for_query : pd.DataFrame,
                            low_concern_descriptions : pd.DataFrame) -> list:
    """
    Performs the logic of a low-concern housekeeping protein dataframe for a
    specific query, and determines whether it contains the criteria to clear 
    a specified hit, for a specified region.

    Benign proteins should have a coverage of at least 80% of the region they are 
    aiming to cover, as well as be at least 50 nucleotide in length.
    """
    logger.debug("\t\t\tChecking query (%s) hit %s for proteins of low-concern",  query.name, hit.name)

    # Ignore this region, if there are no overlapping hits.
    low_concern_protein_for_query_trimmed = _trim_to_region(low_concern_protein_for_query, region).copy()
    if low_concern_protein_for_query_trimmed.empty:
        logger.debug("\t\tNo overlapping Proteins of low-concern for %s", hit.name)
        return []

    low_concern_protein_for_query_trimmed = _calculate_coverage(low_concern_protein_for_query_trimmed, region)
    logger.debug("\tCoverage Data for %s: shape: %s preview:\n%s", 
                     hit.name, low_concern_protein_for_query_trimmed.shape, 
                     low_concern_protein_for_query_trimmed[["coverage_nt", "coverage_ratio"]].head())
    # Ensure that this low_concern hit covers at least 50 nucleotides of the query
    low_concern_protein_for_query_trimmed = low_concern_protein_for_query_trimmed[
        low_concern_protein_for_query_trimmed["coverage_nt"] > MINIMUM_PEPTIDE_COVERAGE]

    if low_concern_protein_for_query_trimmed.empty:
        logger.info("\t --> Insufficient coverage (<%i bp) for housekeeping protein to clear %s (%i-%i).",
        int(MINIMUM_PEPTIDE_COVERAGE), hit.name, region.query_start, region.query_end)
        return []

    # Ensure that this low_concern hit covers at least 80% of the regulated region.
    low_concern_protein_for_query_trimmed = low_concern_protein_for_query_trimmed[
        low_concern_protein_for_query_trimmed["coverage_ratio"] > MINIMUM_QUERY_COVERAGE_FRACTION]

    low_concern_protein_for_query_trimmed = low_concern_protein_for_query_trimmed.reset_index(drop=True)

    if low_concern_protein_for_query_trimmed.empty:
        logger.info("\t --> Insufficient coverage (<%i%%) for housekeeping protein to clear %s (%i-%i).",
        int(MINIMUM_QUERY_COVERAGE_FRACTION*100), hit.name, region.query_start, region.query_end)
        return []

    # Report top hit for Protein / RNA / Synbio
    low_concern_hit = low_concern_protein_for_query_trimmed["subject title"][0]
    low_concern_hit_description = str(*low_concern_descriptions["Description"][low_concern_descriptions["ID"] == low_concern_hit])
    match_ranges = [
        MatchRange(
        float(low_concern_protein_for_query_trimmed['evalue'][0]),
        int(low_concern_protein_for_query_trimmed['s. start'][0]), int(low_concern_protein_for_query_trimmed['s. end'][0]),
        int(low_concern_protein_for_query_trimmed['q. start'][0]), int(low_concern_protein_for_query_trimmed['q. end'][0])
        )
    ]
    low_concern_hit_outcome = HitResult(
            HitScreenStatus(
                ScreenStatus.PASS,
                ScreenStep.LOW_CONCERN_PROTEIN
            ),
            low_concern_hit,
            low_concern_hit_description,
            match_ranges,
            annotations={"Coverage: ":float(low_concern_protein_for_query_trimmed['coverage_ratio'][0])}
        )

    # Rarely, something can be cleared that is already cleared, no need to report on that.
    if hit.recommendation.status not in {ScreenStatus.CLEARED_FLAG, ScreenStatus.CLEARED_WARN}:
        logger.info("\t --> Clearing %s (%s) with house-keeping protein, %s",
                        hit.name, hit.recommendation.status, low_concern_hit_outcome.name)
        hit.recommendation.status = hit.recommendation.status.clear()
        low_concern_hit_outcome.recommendation.status = hit.recommendation.status

    return [low_concern_hit_outcome]

def _filter_low_concern_rna(query : Query,
                            hit : HitResult,
                            region : MatchRange,
                            low_concern_rna_for_query : pd.DataFrame) -> list:
    logger.debug("\t\t\tChecking query (%s) hit %s for RNA of low-concern",  query.name, hit.name)
    
    # Filter low_concern RNA for relevance...
    low_concern_rna_for_query_trimmed = _trim_to_region(low_concern_rna_for_query, region).copy()
    if low_concern_rna_for_query_trimmed.empty:
        logger.debug("\t\t\tNo overlapping low-concern RNA regions for %s (%i-%i)",
                     hit.name, region.query_start, region.query_end)
        return []

    low_concern_rna_for_query_trimmed = _calculate_coverage(low_concern_rna_for_query_trimmed, region)
    logger.debug("\t\t\tCoverage Data for %s: shape: %s preview:\n%s", 
                     hit.name, low_concern_rna_for_query_trimmed.shape, 
                     low_concern_rna_for_query_trimmed[["coverage_nt", "coverage_ratio"]].head())
    
    low_concern_rna_for_query_passed = low_concern_rna_for_query_trimmed[
        (region.length() - low_concern_rna_for_query_trimmed["coverage_nt"]) < MINIMUM_RNA_BASEPAIR_COVERAGE]
    low_concern_rna_for_query_passed = low_concern_rna_for_query_passed.reset_index(drop=True)

    logger.debug("\t\tSummary RNA of low-concern for %s: shape: %s preview:\n%s", 
                hit.name, low_concern_rna_for_query_trimmed.shape,
                low_concern_rna_for_query_trimmed[["coverage_nt", "coverage_ratio"]].head())
    
    if not low_concern_rna_for_query_passed.empty:
        low_concern_hit = low_concern_rna_for_query_trimmed["subject title"][0]
        low_concern_hit_description =  low_concern_rna_for_query_trimmed["description of target"][0]
        match_ranges = [
            MatchRange(
            float(low_concern_rna_for_query_trimmed['evalue'][0]),
            int(low_concern_rna_for_query_trimmed['s. start'][0]), int(low_concern_rna_for_query_trimmed['s. end'][0]),
            int(low_concern_rna_for_query_trimmed['q. start'][0]), int(low_concern_rna_for_query_trimmed['q. end'][0])
            )
        ]
        low_concern_hit_outcome = HitResult(
                HitScreenStatus(
                    ScreenStatus.PASS,
                    ScreenStep.LOW_CONCERN_RNA
                ),
                low_concern_hit,
                low_concern_hit_description,
                match_ranges,
            )
        
        if hit.recommendation.status not in {ScreenStatus.CLEARED_FLAG, ScreenStatus.CLEARED_WARN}:
            logger.info("\t --> Clearing %s %s (region %i-%i), with RNA of low-concern %s",
                        hit.recommendation.status, hit.name, region.query_start, region.query_end, low_concern_hit_outcome.name)
            hit.recommendation.status = hit.recommendation.status.clear()
            low_concern_hit_outcome.recommendation.status = hit.recommendation.status

        return [low_concern_hit_outcome]
    
    logger.info("Clear failed for %s (%s) as RNA of low-concern has >%i bases unaccounted for.",
                    hit.name, hit.recommendation.status, MINIMUM_RNA_BASEPAIR_COVERAGE)
    return []

def _filter_low_concern_dna(query : Query,
                          hit : HitResult,
                          region : MatchRange,
                          low_concern_dna_for_query : pd.DataFrame) -> list:
    logger.debug("\t\t\tChecking query (%s) hit %s for low-concern Synbio",  query.name, hit.name)

    # Filter low_concern SynBio for relevance...
    low_concern_dna_for_query_trimmed = _trim_to_region(low_concern_dna_for_query, region).copy()
    if low_concern_dna_for_query_trimmed.empty:
        logger.debug("\t ... No overlapping low-concern DNA regions for %s (%i-%i)",
                     hit.name, region.query_start, region.query_end)
        return []

    low_concern_dna_for_query_trimmed = _calculate_coverage(low_concern_dna_for_query_trimmed, region)
    logger.debug("\t\t\tCoverage Data for %s: shape: %s preview:\n%s", 
                     hit.name, low_concern_dna_for_query_trimmed.shape, 
                     low_concern_dna_for_query_trimmed[["coverage_nt", "coverage_ratio"]].head())

    low_concern_dna_for_query_trimmed = low_concern_dna_for_query_trimmed[
        low_concern_dna_for_query_trimmed["coverage_ratio"] > MINIMUM_SYNBIO_COVERAGE_FRACTION]
    low_concern_dna_for_query_trimmed = low_concern_dna_for_query_trimmed.reset_index(drop=True)
    logger.debug("\t\t\tFiltered Coverage Synbio for %s: shape: %s preview:\n%s",
                     hit.name, low_concern_dna_for_query_trimmed.shape,
                     low_concern_dna_for_query_trimmed.head())

    if low_concern_dna_for_query_trimmed.empty:
        logger.info("\t --> Synbio sequences <80%% coverage achieved over hit %s over %i-%i.",
                        hit.name, region.query_start, region.query_end)
        return []

    low_concern_hit = low_concern_dna_for_query_trimmed["subject title"][0]
    low_concern_hit_description =  low_concern_dna_for_query_trimmed["subject title"][0]
    match_ranges = [
        MatchRange(
        float(low_concern_dna_for_query_trimmed['evalue'][0]),
        int(low_concern_dna_for_query_trimmed['s. start'][0]), int(low_concern_dna_for_query_trimmed['s. end'][0]),
        int(low_concern_dna_for_query_trimmed['q. start'][0]), int(low_concern_dna_for_query_trimmed['q. end'][0])
        )
    ]
    low_concern_hit_outcome = HitResult(
            HitScreenStatus(
                ScreenStatus.PASS,
                ScreenStep.LOW_CONCERN_DNA
            ),
            low_concern_hit,
            low_concern_hit_description,
            match_ranges,
        )
    
    logger.debug("Processing low-concern Hit: %s", low_concern_hit_outcome)
    
    if hit.recommendation.status not in {ScreenStatus.CLEARED_FLAG, ScreenStatus.CLEARED_WARN}:
        logger.info("\t --> Clearing %s %s region %i-%i as low-concern DNA, with synthetic biology part %s",
                        hit.recommendation.status, hit.name, region.query_start, region.query_end, low_concern_hit_outcome.name)
        hit.recommendation.status = hit.recommendation.status.clear()
        low_concern_hit_outcome.recommendation.status = hit.recommendation.status
    return [low_concern_hit_outcome]

def _update_low_concern_data_for_query(query : Query,
                                  low_concern_protein : pd.DataFrame,
                                  low_concern_rna : pd.DataFrame,
                                  low_concern_dna : pd.DataFrame,
                                  low_concern_descriptions : pd.DataFrame):
    """
    For a single query, look at all three low_concern database outputs, and update the 
    single queries hit descriptions to record all the low_concern hits, as well as clear
    any overlapping WARN or FLAG hits.
    """

    logger.debug("\t\tProcessing low-concern data for query: %s", query.name)

    low_concern_protein_for_query = pd.DataFrame()
    low_concern_rna_for_query = pd.DataFrame()
    low_concern_dna_for_query = pd.DataFrame()

    logger.info(" Low-concern checks for %s", query.name)

    # We only care about the low_concern data for this query, but we need to check for size.
    if not low_concern_protein.empty:
        low_concern_protein_for_query = low_concern_protein[
            low_concern_protein["query name"].str.rsplit("_", n=1).str[0] == query.name
        ]

    if not low_concern_rna.empty:
        low_concern_rna_for_query = low_concern_rna[
            low_concern_rna["query name"] == query.name
        ]

    if not low_concern_dna.empty:
        low_concern_dna_for_query = low_concern_dna[
            low_concern_dna["query acc."] == query.name
        ]

    new_low_concern_hits = []
    # Separated for logging purposes only.
    new_low_concern_protein_hits = []
    new_low_concern_rna_hits = []
    new_low_concern_dna_hits = []

    # Check every region, of every hit that is a FLAG or WARN, against the Benign screen outcomes.
    for hit in query.result.hits.values():
        # Ignore regions that don't require clearing...
        if (hit.recommendation.status not in {
            ScreenStatus.FLAG,
            ScreenStatus.WARN
            }):
            logger.debug("\t\t\tIgnoring %s [%s], nothing to clear.", hit.name, hit.recommendation.status)
            continue
        if (hit.recommendation.from_step == ScreenStep.BIORISK):
            logger.debug("\t\t\tIgnoring %s [%s], cannot clear from Biorisk step.", hit.name, hit.recommendation.status)
            continue

        cleared_regions = 0
        total_regions_to_clear = len(hit.ranges)
        logger.debug("Hit has %i regions required to clear.", total_regions_to_clear)

        for region in hit.ranges:

            if not low_concern_protein_for_query.empty:
                query.mark_as_hit()
                low_concern_proteins = _filter_low_concern_proteins(query, hit, region,
                                            low_concern_protein_for_query,
                                            low_concern_descriptions)
                if low_concern_proteins:
                    cleared_regions += 1
                    new_low_concern_protein_hits.extend(low_concern_proteins)
            
            if not low_concern_rna_for_query.empty:
                query.mark_as_hit()
                low_concern_rna = _filter_low_concern_rna(query, hit, region,
                                       low_concern_rna_for_query)
                if low_concern_rna:
                    cleared_regions += 1
                    new_low_concern_rna_hits.extend(low_concern_rna)
                
            if not low_concern_dna_for_query.empty:
                query.mark_as_hit()
                low_concern_dna = _filter_low_concern_dna(query, hit, region,
                                          low_concern_dna_for_query)
                if low_concern_dna:
                    cleared_regions += 1
                    new_low_concern_dna_hits.extend(low_concern_dna)

        # Quickly check if this hit had all regions cleared.
        if cleared_regions >= total_regions_to_clear:
            logger.debug("Cleared %i / %i regions. This hit is cleared.", cleared_regions, total_regions_to_clear)
            hit.recommendation.status.clear()
        else:
            logger.debug("Only cleared %i/%i regions. This hit was not cleared.", cleared_regions, total_regions_to_clear)
            hit.recommendation.status = hit.recommendation.status.revert_clear()

    # Report Query Specific lack of hits.
    if len(new_low_concern_protein_hits) == 0:
        logger.info("\t ... no housekeeping protein hits.")
    if len(new_low_concern_rna_hits) == 0:
        logger.info("\t ... no low-concern RNA hits.")
    if len(new_low_concern_dna_hits) == 0:
        logger.info("\t ... no low-concern DNA hits.")

    new_low_concern_hits.extend(new_low_concern_protein_hits)
    new_low_concern_hits.extend(new_low_concern_rna_hits)
    new_low_concern_hits.extend(new_low_concern_dna_hits)

    logger.debug("\tNew low-concern hits added: %i", len(new_low_concern_hits))
    for low_concern_addition in new_low_concern_hits:
        query.result.add_new_hit_information(low_concern_addition)
        logger.debug("\t\tAdding low-concern Hit: %s", low_concern_addition)

def parse_low_concern_hits(protein_handler : HmmerHandler,
                            rna_handler : CmscanHandler,
                            dna_handler : BlastNHandler,
                            queries : dict[str,Query],
                            low_concern_desc : pd.DataFrame):
    """
    Parses the outputs from protein (HMMER), RNA (CMSCAN), and DNA (BLASTN) low-concern screens.

    This function processes hit data from the respective handlers, and updates each query's result handle with
    low-concern status information based on these outcomes. 
    Queries with no relevant flagged or warned hits are marked as SKIP. 
    ----
    ## Parameters:
    * `protein_handler` (HmmerHandler): Handler for HMMER-based protein search results.
    * `rna_handler` (CmscanHandler): Handler for CMSCAN-based RNA search results.
    * `dna_handler` (BlastNHandler): Handler for BLASTN-based DNA search results.
    * `queries` (dict[str, Query]): Dictionary of Query objects, i.e. screen fasta input.
    * `low_concern_desc` (pd.DataFrame): DataFrame containing descriptions for low-concern hits.
    -----
    ## Returns:
    * `None`: The function updates `queries` in place.
    """
    # Reading empty outcomes should result in empty DataFrames, not errors.
    low_concern_protein_screen_data = protein_handler.read_output()
    append_nt_querylength_info(low_concern_protein_screen_data, queries)
    recalculate_hmmer_query_coordinates(low_concern_protein_screen_data)
    logger.debug("\tLow-concern Protein Data: shape: %s preview:\n%s",
                 low_concern_protein_screen_data.shape, low_concern_protein_screen_data.head())
    
    low_concern_rna_screen_data = rna_handler.read_output()
    logger.debug("\tLow-concern RNA Data: shape: %s preview:\n%s",
                low_concern_rna_screen_data.shape, low_concern_rna_screen_data.head())
    
    low_concern_dna_screen_data = dna_handler.read_output()
    low_concern_dna_screen_data = get_top_hits(low_concern_dna_screen_data)
    logger.debug("\tLow-concern Synbio Top Hits Data: shape: %s preview:\n%s",
                low_concern_dna_screen_data.shape, low_concern_dna_screen_data.head())
    
    for query in queries.values():

        skip = True
        for hit in query.result.hits.values():
            if hit.recommendation.status in {
                ScreenStatus.FLAG,
                ScreenStatus.WARN
            }:
                skip = False

        if skip:
            query.result.status.low_concern = ScreenStatus.SKIP
        
        if query.result.status.low_concern == ScreenStatus.SKIP:
            logger.debug("Skipping query %s, no regulated regions to clear.", query.name)
            continue
        
        _update_low_concern_data_for_query(query,
                                      low_concern_protein_screen_data,
                                      low_concern_rna_screen_data,
                                      low_concern_dna_screen_data,
                                      low_concern_desc)

        # Calculate the Benign Screen outcomes for each query.
        query.result.status.low_concern = ScreenStatus.PASS
        # If any hits are still warnings, or flags, propagate that to the low_concern step.
        for flagged_hit in query.result.get_flagged_hits():
            query.result.status.update_step_status(
               ScreenStep.LOW_CONCERN_DNA, flagged_hit.recommendation.status
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