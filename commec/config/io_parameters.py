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
import multiprocessing
import importlib.resources

import yaml
from yaml.parser import ParserError

from commec.config.query import Query
from commec.config.constants import DEFAULT_CONFIG_YAML_PATH
from commec.utils.file_utils import expand_and_normalize
from commec.utils.dict_utils import deep_update

logger = logging.getLogger(__name__)

class ScreenIOParameters:
    """
    Container for input settings constructed from arguments to `screen`.
    """
    def __init__(self, args: argparse.Namespace):
        # Inputs that do no have a package-level default, since they are specific to each run
        self.db_dir = args.database_dir
        input_fasta = args.fasta_file
        
        # Set up output files
        self.output_prefix = self.get_output_prefix(input_fasta, args.output_prefix)
        self.output_screen_file = f"{self.output_prefix}.screen"
        self.tmp_log = f"{self.output_prefix}.log.tmp"
        self.output_json = f"{self.output_prefix}.output.json"

        # Query
        self.query: Query = Query(args.fasta_file)

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

        # Sanity check threads settings
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
            self._update_config_from_yaml(cli_config_yaml)

        # Override configuration with any user-provided CLI arguments
        self._update_config_from_cli(args)

        # Update paths in configuration using appropriate string substitution
        self._format_config_paths(self.db_dir)

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
        for arg in args.user_specified_args:
             if arg in self.config and hasattr(args, arg):
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
                    path: '{default}nt_blast/nt'

        This script will update the dictionary to propagate these substitutions.
        If a database directory is provided, it will override the base_path provided in the yaml.
        """
        if self.config.get("base_paths"):
            try:
                base_paths = self.config["base_paths"]
                if db_dir_override is not None:
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
    def get_output_prefix(input_file: str | os.PathLike, prefix_arg="") -> str:
        """
        Returns a prefix that can be used for all output files.

        - If no prefix was given, use the input filename.
        - If a directory was given, use the input filename as file prefix within that directory.
        """
        input_dir = os.path.dirname(input_file)
        # Get file stem (e.g. /home/user/commec/testing_cm_02.fasta -> testing_cm_02)
        input_name = os.path.splitext(os.path.basename(input_file))[0]
        # Take only the first 64 characters of the input name
        if len(input_name) > 64:
            input_name = input_name[:64]

        # If no prefix given, use input filepath without extension
        if not prefix_arg:
            return os.path.splitext(input_file)[0]
        if os.path.isdir(prefix_arg) or prefix_arg.endswith(os.path.sep) or prefix_arg in {".", "..", "~"}:
            # Make the directory if it doesn't exist
            #if not os.path.isdir(prefix_arg):
            os.makedirs(prefix_arg, exist_ok=True)
            # Use the input filename as a prefix within that directory (stripping out the path)
            return os.path.join(prefix_arg, input_name)

        # Existing, non-directory prefixes can be used as-is
        return prefix_arg

    def output_yaml(self, output_filepath : str | os.PathLike):
        """
            Takes the current state of the yaml configuration dictionary and 
            outputs it to a file for posterity, er, reproducibility.
            
            Parameters:
                output_filepath (str | os.PathLike): Path to the output YAML file.
        """
        with open(output_filepath, "w", encoding="utf-8") as stream_out:
            yaml.safe_dump(self.config, stream_out, default_flow_style=False)


    def clean(self):
        """
        Tidy up directories and temporary files after a run.
        """
        if self.config["do_cleanup"]:
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

    @property
    def should_do_protein_screening(self) -> bool:
        return not self.config["in_fast_mode"]

    @property
    def should_do_nucleotide_screening(self) -> bool:
        return not (self.config["in_fast_mode"] or self.config["skip_nt_search"])

    @property
    def should_do_benign_screening(self) -> bool:
        return self.should_do_protein_screening or self.should_do_nucleotide_screening
