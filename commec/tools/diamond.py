#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Handler for Diamond blastx search of protein databases using nucleotide query.
Instantiate a DiamondHnadler, with input local database, input fasta, and output file.
Alter member variables after instantiation to change behaviour, before calling search().
Throws if inputs are invalid. Creates a temporary log file, which is deleted on completion.
"""
import os
import glob
import subprocess
from typing import Optional
import logging
from math import gcd as greatest_common_denominator
from multiprocessing import Pool

from commec.tools.blast_tools import BlastHandler
from commec.tools.search_handler import SearchToolVersion


class DiamondHandler(BlastHandler):
    """
    Diamond blastx search handler. Initialise instance with inputs and outputs,
    alter arguments as desired, and call search() to run Diamond blastx via CLI.
    Concatenates all diamond outputs into a single output file.
    """

    def __init__(self, database_file: str, input_file: str, out_file: str, **kwargs):
        super().__init__(database_file, input_file, out_file, **kwargs)
        self.frameshift: int = 15
        self.do_range_culling = True
        self.jobs: Optional[int] = None
        self.output_format = "6"
        self.output_format_tokens = [
            "qseqid",
            "stitle",
            "sseqid",
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
        ]

        self.db_files = None
        self.find_db_files()

        # Set during `self.search`
        self.concurrent_runs = None
        self.threads_per_run = None

    def find_db_files(self):
        """
        Find all files matching the pattern nr*.dmnd in DB_PATH -- we have found that it improves
        speed to split the nr database into multiple parts and process the chunks in parallel.
        """
        db_suffix = "nr*.dmnd"
        self.db_files = glob.glob(f"{self.db_directory}/{db_suffix}")
        if len(self.db_files) == 0:
            raise FileNotFoundError(
                f"Mandatory Diamond database directory {self.db_directory} "
                f"contains no databases matching the pattern: {db_suffix}"
                " Screening directory can be set via --databases option or --config yaml."
            )

    def run_diamond_search(self, args):
        """Wrapper for run_as_subprocess, used for pooling processes."""
        db_index, db_file = args
        output_file = f"{self.out_file}.{db_index}.tsv"
        log_file = f"{self.temp_log_file}_{db_index}"

        command = [
            "diamond", "blastx", "--quiet",
            "-d", db_file,
            "--threads", str(self.threads_per_run),
            "-q", self.input_file,
            "-o", output_file,
            "--frameshift", str(self.frameshift),
            "--outfmt", self.output_format,
            *self.output_format_tokens
        ]
        if self.do_range_culling:
            command.append("--range-culling")

        self.run_as_subprocess(command, log_file, True)

    def determine_runs_and_threads(
        self, max_threads: int, number_of_databases: int
    ) -> tuple[int, int]:
        """
        Determine the optimal number of Diamond runs and processors per run for best CPU utilization
        and efficiency.
        """
        if self.jobs is not None:
            n_concurrent_runs = self.jobs
        else:
            n_concurrent_runs = greatest_common_denominator(
                number_of_databases, 
                max_threads
                )

        n_threads_per_run = max_threads // n_concurrent_runs

        if n_concurrent_runs < 1:
            logging.info(
                "WARNING: Number of concurrent Diamond runs cannot be < 1. Resetting to 1..."
            )
            n_concurrent_runs = 1

        if n_threads_per_run < 1:
            logging.info(
                "WARNING: Number of threads per Diamond run cannot be < 1. Resetting to 1..."
            )
            n_threads_per_run = 1

        if number_of_databases < n_concurrent_runs:
            logging.info(
                "WARNING: Excessive number of requested concurrent Diamond jobs %i."
                " Resetting to number of Diamond databases %i...",
                n_concurrent_runs,
                number_of_databases,
            )
            n_concurrent_runs = number_of_databases

        return n_concurrent_runs, n_threads_per_run

    def _search(self):
        """
        Search the DIAMOND-formatted nr protein database for matches to the query file.
        """
        self.find_db_files()  # Make sure db files are up to date
        n_diamond_dbs = len(self.db_files)

        self.concurrent_runs, self.threads_per_run = self.determine_runs_and_threads(
            self.threads, n_diamond_dbs
        )
        logging.info(
            "Processing %i Diamond dbs using %i concurrent runs with %i threads per run.",
            n_diamond_dbs,
            self.concurrent_runs,
            self.threads_per_run,
        )

        self.warn_if_nonoptimal_cpu_utilization(n_diamond_dbs)

        # Run each input, with pooling.
        with Pool(self.concurrent_runs) as pool:
            pool.map(self.run_diamond_search, enumerate(self.db_files, 1))

        # Concatenate all output files, and remove the temporary ones.
        temp_outfiles = (f"{self.out_file}.{i}.tsv" for i in range(1, len(self.db_files) + 1))
        self.concatenate_concurrent_outputs(self.out_file, temp_outfiles)
        temp_logfiles = (f"{self.temp_log_file}_{i}" for i in range(1, len(self.db_files) + 1))
        self.concatenate_concurrent_outputs(self.temp_log_file, temp_logfiles)

    def warn_if_nonoptimal_cpu_utilization(self, n_diamond_dbs):
        """
        Let the user know if, based on the number of DIAMOND jobs and threads, the CPU appears
        likely to be over- or under-utlized.
        """
        if self.threads_per_run * self.concurrent_runs < self.threads:
            logging.info(
                "WARNING: Total number of threads across concurrent Diamond job pools [%i*%i]"
                " is less than number of allocated threads [%i]. CPU may be underutilised.",
                self.threads_per_run,
                self.concurrent_runs,
                self.threads,
            )
        if self.threads_per_run * self.concurrent_runs > self.threads:
            logging.info(
                "WARNING: The number of Diamond job pools [%i], each using [%i] threads, may"
                " exceed maximum threads [%i]. CPU may be bottlenecked.",
                self.concurrent_runs,
                self.threads_per_run,
                self.threads,
            )
        if n_diamond_dbs % self.concurrent_runs > 0:
            logging.info(
                "WARNING: The number of Diamond database files [%i] is not "
                " divisible by concurrent jobs [%i]. CPU may be underutilised.",
                n_diamond_dbs,
                self.concurrent_runs,
            )

    def concatenate_concurrent_outputs(self, base_filename, temp_files):
        """
        Parallel Diamond processing produces both logfiles and output files that should be
        concatenated together.
        """
        with open(base_filename, "w", encoding="utf-8") as outfile:
            for temp_file in temp_files:
                if os.path.exists(temp_file):
                    with open(temp_file, "r", encoding="utf-8") as infile:
                        outfile.write(infile.read())
                    os.remove(temp_file)

    def get_version_information(self) -> SearchToolVersion:
        try:
            self.find_db_files()
            result = subprocess.run(
                ["diamond", "dbinfo", "--db", self.db_files[0]],
                capture_output=True,
                text=True,
                check=True,
            )
            lines = result.stdout.splitlines()
            database_info: str = lines[1]

            result = subprocess.run(
                ["diamond", "version"], capture_output=True, text=True, check=True
            )

            tool_info: str = result.stdout.strip()
            return SearchToolVersion(tool_info, database_info)
        except subprocess.CalledProcessError:
            return None
