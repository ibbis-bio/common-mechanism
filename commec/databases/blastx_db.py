#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Defines the `Input and Output Screen Parameters` class, and associated dataclasses.
"""

from commec.databases.database import DatabaseHandler

class BlastXDataBase(DatabaseHandler):
    """ A Database handler specifically for use with BlastX files for commec screening. 
    Allows for full customization of any of the callable blast flags. A better implementation
    may be achieved using import blast"""

    def __init__(self, directory : str, database_file : str, input_file : str, out_file : str):
        super().__init__(directory, database_file, input_file, out_file)
        # We fill this with the defaults, however they can always be overridden after instancing, before screening.
        self.arguments_dictionary = {
            "-num_threads": 1,
            "-evalue": "1e-10",
            "-word_size": 6,
            "-threshold": 21,
            "-max_target_seqs": 5000,
            "-culling_limit": 50,
            "-window_size": 40,
            "-matrix": "BLOSUM62",
            "-gapopen": 11,
            "-gapextend": 1,
            "-seg": "yes",
            "-outfmt": ["7", "qacc",   "stitle",   "sacc",   "staxids",
                        "evalue", "bitscore", "pident", "qlen",
                        "qstart", "qend",     "slen",   "sstart",  "send"]
        }
        self.blastcall = "blastx"

    def screen(self):
        command = [
            self.blastcall,
            "-db", self.db_file,
            "-query", self.input_file,
            "-out", self.out_file,
        ]
        command.extend(self.get_arguments())
        self.run_as_subprocess(command,self.temp_log_file)
