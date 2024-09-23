#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Module for a handlers, specifically for calling Blastn command line interface.
Instantiate a BlastNHandler, with input local database, input fasta, and output file.
Alter instants member variable arguments_dictionary with desired alternatives from the defaults.
Throws if inputs are invalid. Creates a temporary log file, which is deleted on completion.
"""

import subprocess
from commec.tools.blast_tools import BlastHandler
from commec.tools.search_handler import SearchToolVersion

class BlastNHandler(BlastHandler):
    """ 
    A Database handler specifically for use with BlastX files for commec screening. 
    Allows for full customization of any of the callable blast flags.
    """

    def __init__(self, database_file : str, input_file : str, out_file : str):
        super().__init__(database_file, input_file, out_file)
        # We fill this with the defaults, however they can always be overridden after instancing, before screening.
        self.arguments_dictionary = {
            "-outfmt": ["7", "qacc",   "stitle",   "sacc",   "staxids",
                        "evalue", "bitscore", "pident", "qlen",
                        "qstart", "qend",     "slen",   "sstart",  "send"],
                        "-num_threads": 8,
            "-evalue": "10",
            "-max_target_seqs": 50,
            "-culling_limit": 5
        }
        self.blastcall = "blastn"

    def search(self):
        command = [
            self.blastcall,
            "-query", self.input_file,
            "-db", self.db_file,
            "-out", self.out_file
        ]
        command.extend(self.get_arguments())
        self.run_as_subprocess(command,self.temp_log_file)

    def get_version_information(self) -> SearchToolVersion:
        try:
            result = subprocess.run(['blastn', '-version'], capture_output=True, text=True, check=True)
            version_info = " ".join(result.stdout.splitlines())
            return version_info
        except subprocess.CalledProcessError:
            return None
