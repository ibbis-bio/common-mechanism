#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Container class to hold information pertaining to queries from an input fasta file,
as well as derived information, such as translated sequences, whether or not
the query was derived from AA or NT.
"""

import logging
import os
import subprocess

logger = logging.getLogger(__name__)


class Query:
    """
    Query to screen, based on an input FASTA. Self-calculates AA version.
    Future direction to back translate to NT when given AA too.
    """

    def __init__(self, input_fasta_filepath: str):
        self.input_fasta_path = input_fasta_filepath
        self.nt_path: str = ""
        self.aa_path: str = ""

    def setup(self, file_prefix: str):
        """
        Translate or reverse translate query, to be ready in AA or NT format.
        """
        self.nt_path = self.get_cleaned_fasta(file_prefix)
        self.aa_path = f"{file_prefix}.transeq.faa"
        self.translate_query()

    def translate_query(self):
        """Run command transeq, to translate our input sequences."""
        logger.debug("Starting translation of FASTA '%s'. Will output to: %s", self.nt_path, self.aa_path)
        command = ["transeq", self.nt_path, self.aa_path, "-frame", "6", "-clean"]
        result = subprocess.run(command, check=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"Input FASTA '{self.input_fasta_filepath}' could not be translated:\n'{result.stderr}'"
            )

    def get_cleaned_fasta(self, out_prefix):
        """
        Return a FASTA where whitespace (including non-breaking spaces) and illegal characters are
        replaced with underscores.
        """
        cleaned_file = f"{out_prefix}.cleaned.fasta"
        with (
            open(self.input_fasta_path, "r", encoding="utf-8") as fin,
            open(cleaned_file, "w", encoding="utf-8") as fout,
        ):
            for line in fin:
                line = line.strip()
                modified_line = "".join(
                    "_" if c.isspace() or c == "\xc2\xa0" or c == "#" else c for c in line
                )
                fout.write(f"{modified_line}{os.linesep}")
        return cleaned_file
