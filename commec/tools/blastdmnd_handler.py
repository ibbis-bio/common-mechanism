#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Defines the `Input and Output Screen Parameters` class, and associated dataclasses.
"""
import os
import glob
import subprocess

from commec.tools.database_handler import DatabaseHandler, DatabaseVersion

#TODO: Ensure updates from Tessas handling of Diamond are correctly incorperated here.
class DiamondDataBase(DatabaseHandler):
    """ A Database handler specifically for use with Diamond files for commec screening. """

    def __init__(self, directory : str, database_file : str, input_file : str, out_file : str):
        super().__init__(directory, database_file, input_file, out_file)
        self.frameshift : int = 15
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
        logs = []
        output_files = []
        output_lognames = []
        for i, db_file in enumerate(db_files, 1):
            output_file = f"{self.out_file}.{i}.tsv"
            log_i_filename : str = self.temp_log_file + "_" + str(i)
            output_files.append(output_file)
            output_lognames.append(log_i_filename)
            #output_log : str = self.temp_log_file + "_" + str(i)
            command = [
                "diamond", "blastx",
                "--quiet",
                "-d", db_file,
                "--threads", str(self.threads),
                "-q", self.input_file,
                "-o", output_file,
                "--frameshift", str(self.frameshift),
                "--outfmt", self.output_format,
            ]
            command.extend(self.output_format_tokens)

            if self.do_range_culling:
                command.append("--range-culling")

            if self.jobs is not None:
                command.extend(["-j", str(self.jobs)])

            log_i = open(log_i_filename, "w", encoding="utf-8")
            print(" ".join(command))
            process = subprocess.Popen(command, stdout=log_i, stderr=subprocess.STDOUT)
            processes.append(process)
            logs.append(log_i)

        # Wait for all processes to finish
        for process in processes:
            process.wait()

        for log_i in logs:
            log_i.close()

        # Concatenate all output files, and remove the temporary ones.
        with open(self.out_file, 'w', encoding="utf-8") as outfile:
            for o_file in output_files:
                if os.path.exists(o_file):
                    with open(o_file, 'r', encoding="utf-8") as infile:
                        outfile.write(infile.read())
                    os.remove(o_file)

        # Logs are concatenated into a global diamond log.
        with open(self.temp_log_file, 'w', encoding="utf-8") as outlog:
            for log_f in output_lognames:
                if os.path.exists(log_f):
                    with open(log_f, 'r', encoding="utf-8") as infile:
                        outlog.write(infile.read())
                    os.remove(log_f)

