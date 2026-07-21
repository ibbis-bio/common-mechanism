#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Script that concatenates all sequences in a FASTA file.

Usage:
    concat_seqs.py input.fasta
"""

import sys

# read in the file based on the command line argument
filename = sys.argv[1]
# open the file
f = open(filename, "r")
# read the file
lines = f.readlines()
# close the file
f.close()

# use the first line as the sequence ID
seq_id = lines[0].rstrip()

# concatenate all other lines that don't start with '>'
seq = ""
for line in lines[1:]:
    if line.startswith(">"):
        continue
    seq += line.rstrip()

# print the sequence ID and the sequence to a FASTA file
# open output file (filename but with _concat appended before suffix)
out_filename = filename.replace(".fasta", "_concat.fasta")
out_file = open(out_filename, "w")
# write the sequence ID and sequence to the file
out_file.write(">" + seq_id + "\n")
out_file.write(seq + "\n")
