#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Handler for BLASTX search of protein databases using nucleotide queries.
Initialise with local input database, fasta to screen, and output file.
Throws error if inputs are invalid. Creates a temporary log file, which is deleted on completion.
"""

import os
import glob
import subprocess
from commec.tools.blast_tools import BlastHandler
from commec.tools.search_handler import SearchToolVersion, DatabaseValidationError


class BlastXHandler(BlastHandler):
    """
    A search handler specifically for BLASTX command-line during commec screening.
    Modify `arguments_dictionary` to change arguments passed to the CLI.
    """

    def __init__(
        self, database_file: str, input_file: str, out_file: str, **kwargs,
    ):
        super().__init__(database_file, input_file, out_file, **kwargs)
        # We fill this with defaults, however they can always be overridden before screening.
        self.arguments_dictionary = {
            "-task": "blastx-fast",
            "-num_threads": self.threads,
            "-mt_mode": 1,
            "-evalue": 1e-10,
            "-max_target_seqs": 500,
            "-culling_limit": 1,
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
        }
        self.blastcall = "blastx"

    def _validate_db(self):
        """
        BLASTX databases are addressed by a prefix with companion index files:
        single-volume `<prefix>.phr`, multi-volume alias `<prefix>.pal`, or unaliased
        shards `<prefix>.<N>.phr`. Validate the configured prefix points at one of these.
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
                f"No BLASTX database files found for prefix '{self.db_file}'."
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
                ["blastx", "-version"], capture_output=True, text=True, check=True
            )
            tool_info = result.stdout.strip()

            result = subprocess.run(
                ["blastdbcmd", "-info", "-db", self.db_file, "-dbtype", "prot"],
                capture_output=True,
                text=True,
                check=True,
            )
            lines = result.stdout.splitlines()
            database_info: str = lines[5] + lines[3]

            return SearchToolVersion(tool_info, database_info)

        except (subprocess.CalledProcessError, FileNotFoundError):
            return SearchToolVersion()
