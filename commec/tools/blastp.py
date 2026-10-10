#!/usr/bin/env python3
# Copyright (c) 2021-2025 International Biosecurity and Biosafety Initiative for Science
"""
Handler for BLASTP search of protein databases using amino acid queries.
Used in place of BLASTX when the input sequence is already protein.
"""

import glob
import os
import subprocess

from commec.tools.blast_tools import BlastHandler
from commec.tools.search_handler import DatabaseValidationError, SearchToolVersion


class BlastPHandler(BlastHandler):
    """
    A search handler for BLASTP during commec screening of amino acid input.
    Modify arguments_dictionary to change arguments passed to the CLI.
    """

    # Step 3: protein input goes directly to blastp; blastx is skipped

    def __init__(
        self,
        database_file: str,
        input_file: str,
        out_file: str,
        **kwargs,
    ):
        super().__init__(database_file, input_file, out_file, **kwargs)
        self.arguments_dictionary = {
            "-task": "blastp-fast",
            "-num_threads": self.threads,
            "-mt_mode": 1,
            "-evalue": 1e-10,
            "-max_target_seqs": 500,
            "-culling_limit": 1,
            "-outfmt": [
                "6",
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
        }
        self.blastcall = "blastp"

    def _validate_db(self):
        """
        BLASTP databases use the same protein index format as BLASTX.
        Validate that the configured prefix points at a single-volume (.phr),
        multi-volume alias (.pal), or unaliased shards (<prefix>.<N>.phr).
        """
        if not os.path.isdir(self.db_directory):
            raise DatabaseValidationError(
                f"No screening database directory found at: {self.db_directory}."
                " Directory path can be set via --databases option or --config yaml."
            )
        if not (
            os.path.isfile(f"{self.db_file}.phr")
            or os.path.isfile(f"{self.db_file}.pal")
            or glob.glob(f"{self.db_file}.[0-9]*.phr")
        ):
            raise DatabaseValidationError(
                f"No BLASTP database files found for prefix '{self.db_file}'."
                " Expected <prefix>.phr, <prefix>.pal, or <prefix>.<N>.phr in the database"
                " directory. Check the prefix set via --databases or --config yaml matches"
                " the BLAST index files on disk."
            )

    def _search(self):
        command = [
            self.blastcall,
            "-db",
            self.db_file,
            "-query",
            self.input_file,
            "-out",
            self.out_file,
        ]
        command.extend(self.format_args_for_cli())
        self.run_as_subprocess(command, self.temp_log_file)

    def get_version_information(self) -> SearchToolVersion:
        try:
            result = subprocess.run(
                ["blastp", "-version"], capture_output=True, text=True, check=True
            )
            tool_info = result.stdout.strip().replace("\t", " ").replace("\n", " ")

            result = subprocess.run(
                ["blastdbcmd", "-info", "-db", self.db_file, "-dbtype", "prot"],
                capture_output=True,
                text=True,
                check=True,
            )
            lines = result.stdout.splitlines()
            lines = [
                line.strip().replace("\t", " ").replace("\n", " ") for line in lines
            ]
            database_info: str = lines[5] + " " + lines[3]

            return SearchToolVersion(tool_info, database_info)

        except (subprocess.CalledProcessError, FileNotFoundError):
            return SearchToolVersion()
