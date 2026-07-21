#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Split a multi-record FASTA file into files with a set number of sequences per file. Output will be
input.#.fa.

Command-line usage:
    split_fasta.py --i input_fasta -n num_seqs
"""

import os
import sys
import argparse
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", dest="i_file", required=True, help="multi-FASTA file to split"
    )
    parser.add_argument(
        "-n",
        "--num",
        dest="num_seqs",
        type=int,
        required=True,
        help="Number of sequences per file (min)",
    )
    args = parser.parse_args()

    basename = os.path.splitext(args.i_file)[0]
    basename = os.path.basename(basename)
    count_curr = 0
    count_total = 0
    num_splits = 0
    sys.stdout.write("\t%i sequences printed (%i splits)" % (count_total, num_splits))
    sys.stdout.flush()
    for record in SeqIO.parse(args.i_file, "fasta"):
        if count_curr == 0:
            num_splits += 1
            o_file = open(basename + "." + str(num_splits) + ".fa", "w")
        SeqIO.write(record, o_file, "fasta")
        count_curr += 1
        count_total += 1
        if count_total % 10000 == 0:
            sys.stdout.write(
                "\r\t%i sequences printed (%i splits)" % (count_total, num_splits)
            )
            sys.stdout.flush()
        if count_curr == args.num_seqs:
            count_curr = 0
            o_file.close()

    sys.stdout.write("\t%i sequences printed (%i splits)\n" % (count_total, num_splits))
    sys.stdout.flush()


if __name__ == "__main__":
    main()
