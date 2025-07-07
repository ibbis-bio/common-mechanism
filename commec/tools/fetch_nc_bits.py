#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Fetch parts of a query that had no high-quality protein matches for use in nucloetide screening.

Usage:
    fetch_nc_bits.py query_name fasta_file_path
"""
import argparse
import logging
import shutil
import re
import pandas as pd
from Bio import SeqIO
from commec.config.query import Query
from commec.tools.blast_tools import get_high_identity_hits
from commec.tools.search_handler import SearchHandler
from commec.config.result import ScreenStatus

logger = logging.getLogger(__name__)

def _get_ranges_with_no_hits(input_df : pd.DataFrame):
    """
    Get indices not covered by the query start / end ranges in the BLAST results.
    """

    assert "q. start" in input_df.columns, (
        "Column \"q. start\" does not exist for get_ranges_with_no_hits().\n"
        f"Existing columns: {', '.join(input_df.columns)}"
    )

    assert "q. end" in input_df.columns, (
        "Column \"q. end\" does not exist for get_ranges_with_no_hits().\n"
        f"Existing columns: {', '.join(input_df.columns)}"
    )

    assert "query length" in input_df.columns, (
        "Column \"query length\" does not exist for get_ranges_with_no_hits().\n"
        f"Existing columns: {', '.join(input_df.columns)}"
    )

    assert not input_df.empty, "Input dataframe for get_ranges_with_no_hits() is empty."

    unique_hits = input_df.drop_duplicates(subset=["q. start", "q. end"])
    hit_ranges = unique_hits[["q. start", "q. end"]].values.tolist()

    # Sort each pair to ensure that start < end, then sort entire list of ranges
    hit_ranges = sorted([sorted(pair) for pair in hit_ranges])

    nc_ranges : list[tuple[int,int]] = []

    # Include the start if the first hit begins more than 50 bp after the start
    if hit_ranges[0][0] > 50:
        nc_ranges.append((1, hit_ranges[0][0] - 1))

    # Add ranges if there is a noncoding region of >=50 between hits
    for i in range(len(hit_ranges) - 1):
        nc_start = hit_ranges[i][1] + 1  # starts after this hit
        nc_end = hit_ranges[i + 1][0] - 1 # ends before next hit

        if nc_end - nc_start + 1 >= 50:
            nc_ranges.append((nc_start, nc_end))

    # Include the end if the last hit ends more than 50 bp before the end
    query_length = input_df["query length"].iloc[0]
    if query_length - hit_ranges[-1][1] >= 50:
        nc_ranges.append((hit_ranges[-1][1] + 1, int(query_length)))

    return nc_ranges

def _set_no_coding_regions(query : Query):
    """Set the query to be entirely non-coding (i.e. no high-quality protein hits)."""
    query.non_coding_regions.append((1, query.length))

def calculate_noncoding_regions_per_query(
        protein_search_handler : SearchHandler,
        queries : dict[str, Query]
        ):
    """
    Fetch noncoding regions > 50bp for every query, and
    updates the Query dictionary to include non-coding meta-data.
    """
    logger.debug("Checking protein hits in: %s", protein_search_handler.out_file)

    if not protein_search_handler.has_hits():
        logger.info("No protein hits found, screening entire sequence.")
        for query in queries.values():
            _set_no_coding_regions(query)
        return

    protein_hits = get_high_identity_hits(protein_search_handler.out_file)

    query_col = "query acc."

    for query in queries.values():
        protein_hits_for_query = protein_hits[protein_hits[query_col] == query.name].copy()

        if protein_hits_for_query.empty:
            logger.info("No protein hits found for %s, screening entire sequence.", query.name)
            _set_no_coding_regions(query)
            continue

        # Correcting query length in nc coordinate output.
        protein_hits_for_query.loc[:, "q.len"] = query.length

        logger.debug("\t --> Protein hits found for %s, fetching nt regions not covered by a 90%% ID hit or better", query.name)

        ranges_to_screen = _get_ranges_with_no_hits(protein_hits_for_query)
        # if the entire sequence, save regions <50 bases, is covered with protein, skip nt scan
        if not ranges_to_screen:
            logger.info("\t --> no noncoding regions >= 50 bases found for %s, skipping nt scan for query.", query.name)
            query.result.status.nucleotide_taxonomy = ScreenStatus.SKIP
            continue

        # Update the list of start and end non-coding tuples for query.
        query.non_coding_regions.extend(ranges_to_screen)

        ranges_str = ", ".join(f"{start}-{end}" for start, end in ranges_to_screen)
        logger.info("\t --> Identified noncoding regions for query %s: [%s]", query.name, ranges_str)
