#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Defines the `ScreenIO` class and associated dataclasses.
Objects responsible for parsing and interpreting user input for
the screen workflow of commec.
"""

import os
import sys
import glob
import argparse
import logging
import multiprocessing
from pprint import pformat
from pathlib import Path
from Bio import SeqIO
import yaml
from Bio.SeqRecord import SeqRecord

import commec.config.yaml_io as YamlIO
from commec.config.query import Query
from commec.utils.file_utils import expand_and_normalize
from commec.config.constants import (
    MINIMUM_QUERY_LENGTH,
    MAXIMUM_QUERY_LENGTH,
    MAXIMUM_FILENAME_SIZE,
    MAXIMUM_QUERY_NAME_LENGTH,
    VALID_BLAST_MT_MODES,
)

logger = logging.getLogger(__name__)


class ScreenIO:
    """
    Container for input settings constructed from arguments to `screen`.
    """

    def __init__(self, args: argparse.Namespace):
        # Inputs that do no have a package-level default, since they are specific to each run
        self.db_dir = args.database_dir
        self.input_fasta_path = os.path.abspath(args.fasta_file)
        output_prefix = args.output_prefix

        # Output folder hierarchy
        base, outputs, inputs = self._get_output_prefixes(
            self.input_fasta_path, output_prefix
        )
        self.directory_prefix = base
        self.output_prefix = outputs
        self.input_prefix = inputs

        # IO files
        self.output_screen_file = f"{self.directory_prefix}.screen.log"
        self.output_json = f"{self.directory_prefix}.output.json"
        self.nt_path = f"{self.input_prefix}.cleaned.fasta"
        self.aa_path = f"{self.input_prefix}.faa"
        self.nc_path = f"{self.input_prefix}.noncoding.fasta"

        # Get configuration based on defaults and CLI args (including YAML config if supplied)
        self.config = {}
        self._read_config(args)

        # Check whether a .screen output file already exists.
        if os.path.exists(self.output_screen_file) and not (
            self.config["force"] or self.config["resume"]
        ):
            logger.error(
                """Screen output %s already exists.
                Either use a different output location, or use --force or --resume to override.
                Aborting Screen.""",
                self.output_screen_file,
            )
            sys.exit(1)

    def setup(self) -> bool:
        """
        Additional setup once the class has been instantiated (i.e. that requires logs).
        """
        # Make sure the number of threads provided by the user makes sense
        if self.config["threads"] > multiprocessing.cpu_count():
            logger.info(
                "Requested allocated threads [%i] is greater"
                " than the detected CPU count of the hardware[%i].",
                self.config["threads"],
                multiprocessing.cpu_count(),
            )

        if self.config["threads"] < 1:
            raise RuntimeError("Number of allocated threads must be at least 1!")

        if self.config["blast_mt_mode"] not in VALID_BLAST_MT_MODES:
            raise RuntimeError(
                f"blast_mt_mode must be one of {VALID_BLAST_MT_MODES}, "
                f"but got {self.config['blast_mt_mode']!r}"
            )

        # Write a clean FASTA that can be used downstream
        self._write_clean_fasta()

        return True

    def parse_input_fasta(self) -> dict[str, Query]:
        """
        Parse queries from FASTA file.
        """
        records = []
        updated_records = []
        queries = {}

        try:
            with open(self.nt_path, "r", encoding="utf-8") as fasta_file:
                records = list(SeqIO.parse(fasta_file, "fasta"))
        except ValueError as e:
            raise IoValidationError(
                f"Input FASTA file: {self.input_fasta_path} is not a valid fasta file."
            ) from e

        if len(records) == 0:
            raise IoValidationError(
                f"Input FASTA file: {self.input_fasta_path}  contains no records!"
            )

        for record in records:
            try:
                query = Query(record)
                if query.name in queries:
                    raise ValueError(
                        f'Duplicate sequence identifier generated: "{query.name}" from record: {record}\n'
                        f"Ensure that the first {MAXIMUM_QUERY_NAME_LENGTH} characters for each fasta record are unique."
                    )
                queries[query.name] = query
                # Override the original cleaned fasta, with queries above a given length and updated names
                if MINIMUM_QUERY_LENGTH < len(record.seq) <= MAXIMUM_QUERY_LENGTH:
                    # Creating new SeqRecord to avoid overwriting the seq_record object inside query and preserve the original seq id
                    updated_records.append(
                        SeqRecord(record.seq, id=query.name, description="")
                    )
            except Exception as e:
                raise IoValidationError(
                    f"Failed to parse input fasta: {self.nt_path}, {e}"
                ) from e

        with open(self.nt_path, "w", encoding="utf-8") as fasta_file:
            SeqIO.write(updated_records, fasta_file, "fasta")

        return queries

    def clean(self):
        """
        Tidy up directories and temporary files after a run.
        """
        if self.config.do_cleanup:
            for pattern in [
                "*hmmscan",
                "*blastn",
                "faa",
                "*blastx",
                "*.tmp",
            ]:
                for file in glob.glob(f"{self.output_prefix}.{pattern}"):
                    if os.path.isfile(file):
                        os.remove(file)

    def _read_config(self, args: argparse.Namespace):
        """
        Get the configuration for this screen run.

        Configuration is read from multiple sources, which can override each other, according
        to the following hierarchy:

            0. (Lowest-priority) defaults from package-level YAML configuration
            1. Contents of a user-defined YAML file provided using the --config argument
            2. (highest-priority) Configuration provided directly as CLI arguments
        """
        self.config = YamlIO.get_defaults()

        # Import a config yaml file if provided.
        cli_config_yaml = args.config_yaml.strip()
        if cli_config_yaml:
            if not os.path.exists(cli_config_yaml):
                raise FileNotFoundError(f"--config YAML not found: {cli_config_yaml}")
            logger.debug("Overriding defaults in with values from %s", cli_config_yaml)
            self.config = YamlIO.update_config_from_yaml(self.config, cli_config_yaml)

        # Override configuration with any user-provided CLI arguments
        self.config = YamlIO.update_config_from_cli(self.config, args)

        # Override the default base path with database directory from cli.
        base_paths = self.config["base_paths"]
        if self.db_dir is not None:
            base_paths["default"] = Path(expand_and_normalize(self.db_dir)).resolve()
            logger.info(
                "Command line arguments updated base databases directory: %s",
                base_paths["default"],
            )
        else:
            # Otherwise update the default database path if it wasn't defined.
            self.db_dir = base_paths["default"]

        # CLI -d accepts relative paths for shell convenience; normalize to absolute.
        # db_dir_override = expand_and_normalize(self.db_dir) if self.db_dir else None

        # Update paths in configuration using appropriate string substitution
        self.config = YamlIO.format_config_paths(self.config)

        logger.debug("Running Screen with the following parameter set:")
        logger.debug(pformat(self.config))

    @staticmethod
    def _get_output_prefixes(input_file: str | os.PathLike, prefix_arg=None) -> str:
        """
        Returns a set of prefixes that can be used for all output files:
            prefix/name
            prefix/output_name/name
            prefix/input_name/name

        - If no prefix was given, use the input filename as name.
        - If a directory was given, use the input filename
            as file prefix within that directory.
        """
        name = os.path.splitext(os.path.basename(input_file))[0]
        name = name[:MAXIMUM_FILENAME_SIZE]

        directory = prefix_arg or os.path.dirname(input_file)

        # Update the directory/name but only:
        # if the prefix argument was given,
        # if it is not a directory
        if prefix_arg:
            if not (
                os.path.isdir(prefix_arg)
                or prefix_arg.endswith(os.path.sep)
                or prefix_arg in {".", "..", "~"}
            ):
                name = os.path.splitext(os.path.basename(directory))[0]
                directory = os.path.dirname(input_file)

        base = directory
        outputs = os.path.join(directory, f"output_{name}/")
        inputs = os.path.join(directory, f"input_{name}/")

        for path in [base, outputs, inputs]:
            os.makedirs(expand_and_normalize(path), exist_ok=True)

        base_prefix = os.path.join(base, name)
        outputs_prefix = os.path.join(outputs, name)
        inputs_prefix = os.path.join(inputs, name)

        return base_prefix, outputs_prefix, inputs_prefix

    def output_yaml(self, output_filepath: str | os.PathLike):
        """
        Takes the current state of the yaml configuration dictionary and
        outputs it to a file for posterity, er, reproducibility.

        Parameters:
            output_filepath (str | os.PathLike): Path to the output YAML file.
        """
        with open(output_filepath, "w", encoding="utf-8") as stream_out:
            yaml.safe_dump(self.config, stream_out, default_flow_style=False)

    def _write_clean_fasta(self) -> str:
        """
        Write a FASTA in which whitespace (excluding in header), including non-breaking spaces and
        illegal characters are replaced with underscores.
        """

        with (
            open(self.input_fasta_path, "r", encoding="utf-8") as fin,
            open(self.nt_path, "w", encoding="utf-8") as fout,
        ):
            for line in fin:
                line = line.strip()
                if line.startswith(">"):
                    modified_line = "".join(
                        "_" if c == "\xc2\xa0" or c == "#" else c for c in line
                    )
                else:
                    modified_line = "".join(
                        "_" if c.isspace() or c == "\xc2\xa0" or c == "#" else c
                        for c in line
                    )
                fout.write(f"{modified_line}{os.linesep}")

    @property
    def should_do_protein_screening(self) -> bool:
        return not self.config["skip_taxonomy_search"]

    @property
    def should_do_nucleotide_screening(self) -> bool:
        return not (
            self.config["skip_taxonomy_search"] or self.config["skip_nt_search"]
        )

    @property
    def should_do_low_concern_screening(self) -> bool:
        return True


class IoValidationError(ValueError):
    """Custom exception for errors when handling input and output with `ScreenIO`."""
