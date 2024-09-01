#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Defines the `Input and Output Screen Parameters` class, and associated dataclasses.
"""
import os
import glob
import subprocess

from commec.databases.database import DatabaseHandler

#TODO: Ensure updates from Tessas handling of Diamond are correctly incorperated here.
class DiamondDataBase(DatabaseHandler):
    """ A Database handler specifically for use with Diamond files for commec screening. """

    def __init__(self, directory : str, database_file : str, input_file : str, out_file : str):
        super().__init__(directory, database_file, input_file, out_file)
        self.frameshift = 15
        self.do_range_culling = True
        self.threads = 1
        self.jobs = None
        self.output_format = "6"
        self.output_format_tokens = [
            "qseqid",   "stitle",   "sseqid",   "staxids",
            "evalue", "bitscore", "pident", "qlen",
            "qstart", "qend",     "slen",   "sstart",  "send"]

    def get_output_format(self) -> str:
        """ Returns a formatted string version of the output format for blastx"""
        return self.output_format + " ".join(self.output_format_tokens)

    def screen(self):
        # Find all files matching the pattern nr*.dmnd in DB_PATH
        db_files = glob.glob(f"{self.db_directory}/nr*.dmnd")

        if len(db_files) == 0:
            raise FileNotFoundError(f"Mandatory Diamond database directory {self.db_directory} contains no databases!")

        # Run diamond blastx in parallel using subprocess
        processes = []
        for i, db_file in enumerate(db_files, 1):
            output_file = f"{self.out_file}.{i}.tsv"
            command = [
                "diamond", "blastx",
                "--quiet",
                "-d", db_file,
                "--threads", str(self.threads),
                "-q", self.input_file,
                "-o", output_file,
                "--outfmt", self.get_output_format(),
                "--frameshift", str(self.frameshift),
            ]

            if self.do_range_culling:
                command.append("--range-culling")

            if self.jobs is not None:
                command.extend(["-j", str(self.jobs)])

            process = subprocess.Popen(command, self.temp_log_file)
            processes.append(process)

        # Wait for all processes to finish
        for process in processes:
            process.wait()

        # Concatenate all output files
        with open("self.out_file", 'w', encoding="utf-8") as outfile:
            for i in range(1, len(db_files) + 1):
                with open(f"{self.out_file}.{i}.tsv", 'r', encoding="utf-8") as infile:
                    outfile.write(infile.read())

        # Remove the individual output files
        for i in range(1, len(db_files) + 1):
            os.remove(f"{self.out_file}.{i}.tsv")
