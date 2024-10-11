#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Handler for BLASTN search of nucleotide databases using nucleotide queries.
Initialise with local input database, fasta to screen, and output file.
Throws error if inputs are invalid. Creates a temporary log file, which is deleted on completion.
"""

import subprocess
from commec.tools.blast_tools import BlastHandler
from commec.tools.search_handler import SearchToolVersion


class BlastNHandler(BlastHandler):
    """
    A search handler specifically for BLASTN command-line during commec screening.
    Modify `arguments_dictionary` to change passed to the command line call.
    """

    def __init__(self, database_file: str, input_file: str, out_file: str):
        super().__init__(database_file, input_file, out_file)
        # We fill this with defaults, however they can always be overridden before screening.
        self.arguments_dictionary = {
            "-outfmt": [
                "7",
                "qacc",
                "stitle",
                "sacc",
                "staxids",
                "evalue",
                "bitscore",
                "pident",
                "qlen",
                "qstart",
                "qend",
                "slen",
                "sstart",
                "send",
            ],
            "-num_threads": 8,
            "-evalue": "10",
            "-max_target_seqs": 50,
            "-culling_limit": 5,
        }
        self.blastcall = "blastn"

    def search(self):
        command = [
            self.blastcall,
            "-query",
            self.input_file,
            "-db",
            self.db_file,
            "-out",
            self.out_file,
        ]
        command.extend(self.format_args_for_cli())
        self.run_as_subprocess(command, self.temp_log_file)

    def get_version_information(self) -> SearchToolVersion:
        try:
            result = subprocess.run(
                ["blastn", "-version"], capture_output=True, text=True, check=True
            )
            version_info = result.stdout.strip()
            return version_info
        except subprocess.CalledProcessError:
            return None
