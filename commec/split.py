#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Split a multi-record FASTA file into individual files, one for each record.

Command-line usage:
    split.py input.fasta
"""

import argparse
import os
import string

from Bio import SeqIO

from commec.utils.file_utils import file_arg

VALID_FILENAME_CHARS = f"-._{string.ascii_letters}{string.digits}"
DESCRIPTION = (
    "Split a multi-record FASTA file into individual files, one for each record"
)


def add_args(parser):
    """
    Add module arguments to an ArgumentParser object.
    """
    parser.add_argument(
        action="store", dest="fasta_file", type=file_arg, help="Input fasta file"
    )
    return parser


def clean_description(description):
    """
    Cleans the description from a sequence record for use as part of a filename.
    """
    cleaned = description.strip()
    cleaned = "".join(x for x in cleaned if x in VALID_FILENAME_CHARS)
    if len(cleaned) > 150:
        cleaned = cleaned[:150]
    return cleaned


def write_split_fasta(fasta_file):
    """
    Parse all sequence records in an input FASTA file, and write a new file for each record.
    """
    output_dir = os.path.dirname(fasta_file)
    fasta_name = os.path.splitext(os.path.basename(fasta_file))[0]

    with open(fasta_file, "r", encoding="utf-8") as input_file:
        for i, record in enumerate(SeqIO.parse(input_file, "fasta")):
            desc = clean_description(record.description)

            # Handle empty descriptions and avoid overwriting input files
            if not desc or desc == fasta_name:
                output_basename = f"{fasta_name}-split-{i}.fasta"
            else:
                output_basename = f"{desc}.fasta"

            output_path = os.path.join(output_dir, output_basename)
            with open(output_path, "w", encoding="utf-8") as output_file:
                output_file.write(f">{desc}{os.linesep}")
                output_file.write(str(record.seq))


def run(parsed_args):
    """
    Wrapper so that args be parsed in main() or commec.py interface.
    """
    write_split_fasta(parsed_args.fasta_file)


def main():
    """
    Main function. Passes FASTA file to `write_split_fasta`.

    Arguments:
      - fasta_file: Path to the input FASTA file.

    """
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
