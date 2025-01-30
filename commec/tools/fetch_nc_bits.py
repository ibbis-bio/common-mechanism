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
from Bio import SeqIO
from commec.tools.blast_tools import get_high_identity_matches
from commec.tools.search_handler import SearchHandler


def get_ranges_with_no_hits(blast_df):
    """
    Get indices not covered by the query start / end ranges in the BLAST results.
    """
    unique_hits = blast_df.drop_duplicates(subset=["q. start", "q. end"])
    hit_ranges = unique_hits[["q. start", "q. end"]].values.tolist()

    # Sort each pair to ensure that start < end, then sort entire list of ranges
    hit_ranges = sorted([sorted(pair) for pair in hit_ranges])

    nc_ranges = []

    # Include the start if the first hit begins more than 50 bp after the start
    if hit_ranges[0][0] > 50:
        nc_ranges.append([1, hit_ranges[0][0] - 1])

    # Add ranges if there is a noncoding region of >=50 between hits
    for i in range(len(hit_ranges) - 1):
        nc_start = hit_ranges[i][1] + 1  # starts after this hit
        nc_end = hit_ranges[i + 1][0] - 1 # ends before next hit

        if nc_end - nc_start + 1 >= 50:
            nc_ranges.append([nc_start, nc_end])

    # Include the end if the last hit ends more than 50 bp before the end
    query_length = blast_df["query length"][0]
    if query_length - hit_ranges[-1][1] >= 50:
        nc_ranges.append([hit_ranges[-1][1] + 1, int(query_length)])

    return nc_ranges


def get_records(fasta_file_path):
    """Parse SeqIO records from input FASTA."""
    with open(fasta_file_path, "r", encoding="utf-8") as fasta_file:
        records = list(SeqIO.parse(fasta_file, "fasta"))
        return records


def write_nc_sequences(nc_ranges, record, outfile):
    """
    Write a FASTA containing only the parts of the record without a high-quality protein match.
    """
    nc_sequences = []
    for start, stop in nc_ranges:
        # subtract 1 just from `start` to adjust to 0-based index and capture full range
        sequence = record.seq[int(start) - 1 : int(stop)]
        # when parsed from a FASTA file, description includes the id:
        # https://biopython.org/docs/latest/Tutorial/chapter_seq_annot.html#seqrecord-objects-from-fasta-files
        nc_sequences.append(f">{record.description} {start}-{stop}\n{sequence}\n")

    with open(outfile, "w", encoding="utf-8") as output_file:
        output_file.writelines(nc_sequences)


def fetch_noncoding_regions(protein_results, query_fasta):
    """Fetch noncoding regions > 50bp and write to a new file."""
    outfile = re.sub(".nr.*", "", protein_results) + ".noncoding.fasta"

    logging.info("Checking protein hits in: %s", protein_results)

    if SearchHandler.is_empty(protein_results) or not SearchHandler.has_hits(
        protein_results
    ):
        logging.info("\t...no protein hits found, screening entire sequence\n")
        shutil.copyfile(query_fasta, outfile)
        return

    blast_df = get_high_identity_matches(protein_results)

    if blast_df.empty:
        logging.info(
            "Protein hits all low percent identity (<90%%) - screening entire sequence"
        )
        shutil.copyfile(query_fasta, outfile)
        return

    query_col = "query acc."
    if blast_df[query_col].nunique() > 1:
        first_query = blast_df[query_col].iloc[0]
        logging.info(
            "WARNING: Only fetching nucleotides from first query [%s] in multi-query results: %s",
            first_query,
            protein_results,
        )
        blast_df = blast_df[blast_df[query_col] == first_query]

    logging.info(
        "Protein hits found, fetching nt regions not covered by a 90%% ID hit or better"
    )
    ranges_to_screen = get_ranges_with_no_hits(blast_df)

    # if the entire sequence, save regions <50 bases, is covered with protein, skip nt scan
    if not ranges_to_screen:
        logging.info(
            "\t\t --> no noncoding regions >= 50 bases found, skipping nt scan\n"
        )
        return

    records = get_records(query_fasta)

    if len(records) > 1:
        logging.info(
            "WARNING: Only fetching nucleotides from first record in multifasta: %s",
            query_fasta,
        )

    ranges_str = ", ".join(f"{start}-{end}" for start, end in ranges_to_screen)
    logging.info("Writing noncoding regions [%s] to: %s", ranges_str, outfile)
    write_nc_sequences(ranges_to_screen, records[0], outfile)


def main():
    """
    Wrapper for parsing arguments direction to fetch_nc_bits if called as main.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(dest="protein_results", help="Results of a protein search")
    parser.add_argument(
        dest="fasta_file_path",
        help="FASTA file from which to fetch regions with no protein hits",
    )
    args = parser.parse_args()
    fetch_noncoding_regions(args.protein_results, args.fasta_file_path)


if __name__ == "__main__":
    main()