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
import importlib.resources
from pprint import pformat

import yaml
from yaml.parser import ParserError

from Bio import SeqIO

from commec.config.query import Query
from commec.config.constants import DEFAULT_CONFIG_YAML_PATH, MINIMUM_QUERY_LENGTH
from commec.utils.file_utils import expand_and_normalize
from commec.utils.dict_utils import deep_update

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
        base, outputs, inputs = self._get_output_prefixes(self.input_fasta_path, output_prefix)
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
            self.config["force"] or self.config["resume"]):
            logger.warning(
                f"""Screen output {self.output_screen_file} already exists.
                Either use a different output location, or use --force or --resume to override.
                Aborting Screen."""
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

        if (
            self.config["diamond_jobs"] is not None
            and self.config["protein_search_tool"] == "blastx"
        ):
            logger.warning(
                "--jobs is a diamond only parameter! Specifying -j (--jobs) without also"
                " specifying -p (--protein-search-tool) as 'diamond' will have no effect!"
            )

        # Write a clean FASTA that can be used downstream
        self._write_clean_fasta()
        return True


    def parse_input_fasta(self) -> dict[str, Query]:
        """
        Parse queries from FASTA file.
        """
        records = []
        queries = {}

        try:
            with open(self.nt_path, "r", encoding = "utf-8") as fasta_file:
                records = list(SeqIO.parse(fasta_file, "fasta"))
        except ValueError as e:
            raise IoValidationError(f"Input FASTA file: {self.input_fasta_path} "
                                    "is not a valid fasta file.") from e

        if len(records) == 0:
            raise IoValidationError(f"Input FASTA file: {self.input_fasta_path} "
                                    " contains no records!")

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
            
        # Don't write a cleaned fasta for queries below a given length.
        records = [record for record in records if len(record.seq) > MINIMUM_QUERY_LENGTH]

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

    def _read_config(self, args: argparse.Namespace):
        """
        Get the configuration for this screen run.

        Configuration is read from multiple sources, which can override each other, according
        to the following hierarchy:

            0. (Lowest-priority) defaults from package-level YAML configuration
            1. Contents of a user-defined YAML file provided using the --config argument
            2. (highest-priority) Configuration provided directly as CLI arguments
        """
        # Read package-level configuration defaults
        default_yaml = importlib.resources.files("commec").joinpath(DEFAULT_CONFIG_YAML_PATH)
        if default_yaml.exists():
            self.config = self._load_config_from_yaml(str(default_yaml))
        else:
            raise FileNotFoundError(
                f"No default yaml found. Expected at {DEFAULT_CONFIG_YAML_PATH}"
                )

        # Override configuration with any in user-provided YAML file
        cli_config_yaml=args.config_yaml.strip()
        if os.path.exists(cli_config_yaml):
            logger.debug(f"Overriding defaults in {default_yaml} with values from {cli_config_yaml}")
            self._update_config_from_yaml(cli_config_yaml)

        # Override configuration with any user-provided CLI arguments
        self._update_config_from_cli(args)

        # Update paths in configuration using appropriate string substitution
        self._format_config_paths(self.db_dir)

        logger.debug("Running Screen with the following parameter set:")
        logger.debug(pformat(self.config))

    def _load_config_from_yaml(self, config_filepath: str | os.PathLike) -> dict:
        """
        Loads a yaml file, ensuring it's a dictionary.
        """
        try:
            with open(config_filepath, "r", encoding = "utf-8") as file:
                config_from_yaml = yaml.safe_load(file)
        except ParserError:
            raise ValueError(f"Invalid yaml syntax in configuration file: {config_filepath}")

        if not isinstance(config_from_yaml, dict):
            raise TypeError(f"Loaded configuration file did not result in a dictionary: {file}")
        return config_from_yaml

    def _update_config_from_yaml(self, config_filepath: str | os.PathLike) -> None:
        """
        Override YAML configuration based on provided YAML file. Items in the provided file, but
        not in the default YAML, will be ignored.
        """
        config_from_yaml = self._load_config_from_yaml(config_filepath)
        self.config = deep_update(self.config, config_from_yaml)

    def _update_config_from_cli(self, args: argparse.Namespace):
        """ 
        Override YAML configuration based on arguments given in the command line.
        Need to reference `user_specified_args` because CLI defaults should not override YAML.
        """
        if not hasattr(args, "user_specified_args"):
            raise ValueError(
                "Missing required 'user_specified_args' in arguments namespace. "
            )

        # Update the YAML default values in the configuration dictionary
        logger.debug("Using the following CLI configuration arguments:")
        logger.debug(pformat(args.user_specified_args))

        for arg in args.user_specified_args:
            if arg in self.config and hasattr(args, arg):
                logger.debug(f"Command line arguments updated '{arg}' to: {getattr(args,arg)}")
                self.config[arg] = getattr(args, arg)

    def _format_config_paths(self,
        db_dir_override: str | os.PathLike = None
    ):
        """
        The YAML file is expected to contain a 'base_paths' key that is referenced in string
        substitutions, so that base paths do not need to be defined more than once. For example:

            base_paths:
                default: path/to/databases/
            databases:
                regulated_nt:
                    path: '{default}nt_blast/core_nt'

        This script will update the dictionary to propagate these substitutions.
        If a database directory is provided, it will override the base_path provided in the yaml.
        """
        if self.config.get("base_paths"):
            try:
                base_paths = self.config["base_paths"]
                if db_dir_override is not None:
                    logger.debug(f"Command line arguments updated base databases directory: {db_dir_override}")
                    base_paths["default"] = db_dir_override
                else:
                    self.db_dir = base_paths["default"]

                # Ensure all the base paths end with a separator
                for key, value in base_paths.items():
                    base_paths[key] = os.path.join(value,'')

                def recursive_format(nested_yaml, base_paths):
                    """
                    Recursively apply string formatting to read paths from nested yaml config dicts.
                    """
                    if isinstance(nested_yaml, dict):
                        return {key : recursive_format(value, base_paths) 
                                for key, value in nested_yaml.items()}
                    if isinstance(nested_yaml, str):
                        try:
                            return nested_yaml.format(**base_paths)
                        except KeyError as e:
                            raise ValueError(
                                f"Unknown base path key referenced in path: {nested_yaml}"
                            ) from e
                    return nested_yaml

                self.config = recursive_format(self.config, base_paths)
            except TypeError:
                pass

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
        name = name[:64]

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

        base_prefix = os.path.join(base,name)
        outputs_prefix =  os.path.join(outputs,name)
        inputs_prefix =  os.path.join(inputs,name)

        return base_prefix, outputs_prefix, inputs_prefix

    def output_yaml(self, output_filepath : str | os.PathLike):
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
        return not self.config["in_fast_mode"]

    @property
    def should_do_nucleotide_screening(self) -> bool:
        return not (self.config["in_fast_mode"] or self.config["skip_nt_search"])

    @property
    def should_do_benign_screening(self) -> bool:
        return True


class IoValidationError(ValueError):
    """Custom exception for errors when handling input and output with `ScreenIO`."""