#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Defines the `ScreenIOParameters` class and associated dataclasses.
Objects responsible for parsing and interpreting user input for
the screen workflow of commec.
"""
import os
import sys
import glob
import argparse
import logging
from dataclasses import dataclass
from typing import Optional
import multiprocessing

import yaml
from yaml.parser import ParserError

from Bio import SeqIO

from commec.config.query import Query
from commec.config.constants import DEFAULT_CONFIG_YAML_PATH

@dataclass
class ScreenConfig:
    """
    Namespace for optional input parameters for screening; provided by parsing command line
    arguments. Default values, where applicable, are stored here.
    """
    threads: int = 1
    protein_search_tool: str = "blastx"
    in_fast_mode: bool = False
    skip_nt_search: bool = False
    do_cleanup: bool = False
    diamond_jobs: Optional[int] = None
    config_yaml_file: str | os.PathLike = DEFAULT_CONFIG_YAML_PATH
    force: bool = False
    resume: bool = False
    benchmark : bool = False

class ScreenIO:
    """
    Container for input settings constructed from arguments to `screen`.
    """

    def __init__(self, args: argparse.ArgumentParser):
        # Process command-line inputs from user
        self.config: ScreenConfig = ScreenConfig(
            threads=args.threads,
            protein_search_tool=args.protein_search_tool,
            in_fast_mode=args.fast_mode,
            skip_nt_search=args.skip_nt_search,
            do_cleanup=args.cleanup,
            diamond_jobs=args.diamond_jobs,
            config_yaml_file=args.config_yaml.strip(),
            force=args.force,
            resume=args.resume,
            benchmark=args.benchmark,
        )

        self.input_fasta_path = args.fasta_file

        # Set up output files
        self.output_prefix = self._get_output_prefix(self.input_fasta_path, args.output_prefix)
        self.nt_path = f"{self.output_prefix}.cleaned.fasta"
        self.aa_path = f"{self.output_prefix}.faa"
        self.output_screen_file = f"{self.output_prefix}.screen"
        self.tmp_log = f"{self.output_prefix}.log.tmp"
        self.output_json = f"{self.output_prefix}.output.json"

        # Parse paths to input databases (may be from CLI or YAML configuration file)
        self.db_dir = args.database_dir
        self.yaml_configuration = {}

        if os.path.exists(self.config.config_yaml_file):
            self._get_configurations_from_yaml(self.config.config_yaml_file, self.db_dir)
        else:
            raise FileNotFoundError(
                "No configuration yaml found. If using a custom file, check the path is correct: "
                + self.config.config_yaml_file
            )

        # Check whether a .screen output file already exists.
        if os.path.exists(self.output_screen_file) and not (
            self.config.force or self.config.resume
        ):
            # Print statement must be used as logging not yet instantiated
            print(
                f"Screen output {self.output_screen_file} already exists. \n"
                "Either use a different output location, or use --force or --resume to override. "
                "\nAborting Screen."
            )
            sys.exit(1)

    def setup(self) -> bool:
        """
        Additional setup once the class has been instantiated (i.e. that requires logs).
        """
        # Make sure the number of threads provided by the user makes sense
        if self.config.threads > multiprocessing.cpu_count():
            logging.info(
                "Requested allocated threads [%i] is greater"
                " than the detected CPU count of the hardware[%i].",
                self.config.threads,
                multiprocessing.cpu_count(),
            )
        if self.config.threads < 1:
            raise RuntimeError("Number of allocated threads must be at least 1!")

        if (
            self.config.diamond_jobs is not None
            and self.config.protein_search_tool == "blastx"
        ):
            logging.info(
                "WARNING: --jobs is a diamond only parameter! "
                "Specifying -j (--jobs) without also specifying "
                "-p (--protein-search-tool), the protein search "
                'tool as "diamond" will have no effect!'
            )

        # Write a clean FASTA that can be used downstream
        self._write_clean_fasta()
        return True

    def parse_input_fasta(self) -> dict[str, Query]:
        """
        Parse queries from FASTA file.
        """
        records = []

        with open(self.nt_path, "r", encoding = "utf-8") as fasta_file:
            queries = {}
            records = list(SeqIO.parse(fasta_file, "fasta"))

        for record in records:
            try:
                query = Query(record)
                if query.name in queries:
                    raise ValueError(f"Duplicate sequence identifier found: {query.name}")
                queries[query.name] = query
                # Override the original cleaned fasta, with updated names.
                record.id = query.name
                record.name = ""
                record.description = ""
            except Exception as e:
                raise IoValidationError(f"Failed to parse input fasta: {self.nt_path}") from e

        with open(self.nt_path, "w", encoding = "utf-8") as fasta_file:
            SeqIO.write(records, fasta_file, "fasta")

        return queries

    def clean(self):
        """
        Tidy up directories and temporary files after a run.
        """
        if self.config.do_cleanup:
            for pattern in [
                "reg_path_coords.csv",
                "*hmmscan",
                "*blastn",
                "faa",
                "*blastx",
                "*dmnd",
                "*.tmp",
            ]:
                for file in glob.glob(f"{self.output_prefix}.{pattern}"):
                    if os.path.isfile(file):
                        os.remove(file)

    def _get_configurations_from_yaml(
        self, config_filepath: str | os.PathLike, db_dir_override: str | os.PathLike = None
    ):
        """
        Read the contents of YAML file configuration.

        The YAML file is expected to contain a 'base_paths' key that is referenced in string
        substitutions, so that base paths do not need to be defined more than once. For example:

            base_paths:
                default: path/to/databases/
            databases:
                regulated_nt:
                    path: '{default}nt_blast/nt'
        """
        try:
            with open(config_filepath, "r", encoding="utf-8") as file:
                config_from_yaml = yaml.safe_load(file)
        except ParserError as e:
            raise ValueError(
                f"Invalid YAML syntax in configuration file: {config_filepath}"
            ) from e

        # Extract base paths for substitution
        missing_sections = {"base_paths", "databases"} - set(config_from_yaml.keys())
        if missing_sections:
            raise ValueError(
                f"Configuration missing required sections: {missing_sections}"
            )

        try:
            base_paths = config_from_yaml["base_paths"]
            if db_dir_override is not None:
                base_paths["default"] = db_dir_override
            else:
                self.db_dir = base_paths["default"]

            def recursive_format(d, base_paths):
                """
                Recursively apply string formatting to read paths from nested yaml config dicts.
                """
                if isinstance(d, dict):
                    return {k: recursive_format(v, base_paths) for k, v in d.items()}
                if isinstance(d, str):
                    try:
                        return d.format(**base_paths)
                    except KeyError as e:
                        raise ValueError(
                            f"Unknown base path key referenced in path: {d}"
                        ) from e
                return d

            config_from_yaml = recursive_format(config_from_yaml, base_paths)
        except TypeError:
            pass

        self.yaml_configuration = config_from_yaml

    @staticmethod
    def _get_output_prefix(input_file, prefix_arg=""):
        """
        Returns a prefix that can be used for all output files.
        - If no prefix was given, use the input filename.
        - If a directory was given, use the input filename 
            as file prefix within that directory.
        """
        if not prefix_arg:
            return os.path.splitext(input_file)[0]
        if (
            os.path.isdir(prefix_arg)
            or prefix_arg.endswith(os.path.sep)
            or prefix_arg in {".", "..", "~"}
        ):
            os.makedirs(prefix_arg, exist_ok=True)
            # Use the input filename as a prefix within that directory (stripping out the path)
            input_name = os.path.splitext(os.path.basename(input_file))[0]
            return prefix_arg + input_name
        # Existing, non-path prefixes can be used as-is
        return prefix_arg

    def _write_clean_fasta(self) -> str:
        """
        Write a FASTA in which whitespace (including non-breaking spaces) and 
        illegal characters are replaced with underscores.
        """
        with (
            open(self.input_fasta_path, "r", encoding="utf-8") as fin,
            open(self.nt_path, "w", encoding="utf-8") as fout,
        ):
            for line in fin:
                line = line.strip()
                modified_line = "".join(
                    "_" if c.isspace() or c == "\xc2\xa0" or c == "#" else c
                    for c in line
                )
                fout.write(f"{modified_line}{os.linesep}")

    @property
    def should_do_protein_screening(self) -> bool:
        return not self.config.in_fast_mode

    @property
    def should_do_nucleotide_screening(self) -> bool:
        return not (self.config.in_fast_mode or self.config.skip_nt_search)

    @property
    def should_do_benign_screening(self) -> bool:
        return True


class IoValidationError(ValueError):
    """Custom exception for errors when handling input and output with `ScreenIO`."""
