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
import re
import argparse
import logging
import multiprocessing
import importlib.resources

import yaml
from yaml.parser import ParserError

from commec.config.query import Query
from commec.config.constants import DEFAULT_CONFIG_YAML_PATH

class ScreenIOParameters:
    """
    Container for input settings constructed from arguments to `screen`.
    """
    def __init__(self, args: argparse.ArgumentParser):
        # Non-yaml Inputs
        cli_config_yaml_filepath=args.config_yaml.strip()
        self.db_dir = args.database_dir
        input_fasta = args.fasta_file
        output_prefix = args.output_prefix

        # Output folder hierarchy
        base, outputs, inputs = self._get_output_prefixes(input_fasta, output_prefix)
        self.directory_prefix = base
        self.output_prefix = outputs
        self.input_prefix = inputs

        self.output_screen_file = f"{self.directory_prefix}.screen"
        self.tmp_log = f"{self.directory_prefix}.log.tmp"

        # Query
        self.query: Query = Query(args.fasta_file)

        # Configuration initialisation:
        self.config = {}

        # Initialize config with package defaults:
        default_config_resource = importlib.resources.files("commec").joinpath(DEFAULT_CONFIG_YAML_PATH)
        if default_config_resource.exists():
            self._lazy_update_from_yaml(str(default_config_resource))
        else:
            raise FileNotFoundError(
                f"No default yaml found. Expected at {DEFAULT_CONFIG_YAML_PATH}"
                )
            
        # Override configuration with provided yaml:
        if os.path.exists(cli_config_yaml_filepath):
            self._lazy_update_from_yaml(cli_config_yaml_filepath)

        # TODO: Consider looking for additional config yaml if none provided at /database/
        
        # Override config with CLI.
        # (If the database directory is provided via CLI, 
        #  use it to override base paths.)
        self._update_config_from_cli(args)

        # Propagate base paths of the configuration:
        self._finalize_config(self.db_dir)

        # Check whether a .screen output file already exists.
        if os.path.exists(self.output_screen_file) and not (
            self.config["force"] or self.config["resume"]):
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
        # Sanity checks on thread input.
        if self.config["threads"] > multiprocessing.cpu_count():
            logging.info(
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
            logging.info(
                "WARNING: --jobs is a diamond only parameter! "
                "Specifying -j (--jobs) without also specifying "
                "-p (--protein-search-tool), the protein search "
                'tool as "diamond" will have no effect!'
            )

        self.query.setup(self.input_prefix)
        return True


    def _update_config_from_cli(self, args: argparse.ArgumentParser):
        """ Maps any CLI, that can update the yaml configuration dictionary, and does so."""
        # Filter to include only explicitly provided arguments
        explicit_args = {k: v for k, v in vars(args).items() if f"--{k.replace('_', '-')}" in sys.argv}

        # Map argparse keys to YAML config keys, some are the same, and can be ignored.
        new_key_names = {
            "fast_mode" : "in_fast_mode",
            "cleanup" : "do_cleanup",
        }
        # Rename arguments based on the Config.yaml mapping
        renamed_args = {new_key_names.get(k, k): v for k, v in explicit_args.items()}

        # Update the configuration dictionary
        for key, value in renamed_args.items():
            self.config[key] = value


    def _lazy_update_from_yaml(self, config_filepath: str | os.PathLike) -> None:
        """
        Loads a yaml file as a dictionary into the Config dictionary.
        The load is done lazily - and basepaths are not parsed and replaced throughout.
        """        
        raw_config_from_yaml = {}
        try:
            with open(config_filepath, "r", encoding = "utf-8") as file:
                raw_config_from_yaml = yaml.safe_load(file)
        except ParserError:
            raise ValueError(f"Invalid yaml syntax in configuration file: {config_filepath}")
        assert isinstance(raw_config_from_yaml, dict), "Loaded configuration file (yaml) did not result in a Dictionary!"
        self.config.update(raw_config_from_yaml)


    def _finalize_config(self,
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

                def recursive_format(dictionary, base_paths):
                    """
                    Recursively apply string formatting to read paths from nested yaml config dicts.
                    """
                    if isinstance(dictionary, dict):
                        return {key : recursive_format(value, base_paths) 
                                for key, value in dictionary.items()}
                    if isinstance(dictionary, str):
                        try:
                            return dictionary.format(**base_paths)
                        except KeyError as e:
                            raise ValueError(
                                f"Unknown base path key referenced in path: {dictionary}"
                            ) from e
                    return dictionary

                self.config = recursive_format(self.config, base_paths)
            except TypeError:
                pass

    def _get_output_prefixes(self, input_file, prefix_arg=""):
        """
        Returns a set of prefixes that can be used for all output files:
            prefix/name
            prefix/output_name/name
            prefix/input_name/name

        - If no prefix was given, use the input filename.
        - If a directory was given, use the input filename as file prefix within that directory.
        """
        name = os.path.splitext(os.path.basename(input_file))[0]
        name = self._decimate_name(name)

        prefix = prefix_arg or name

        base = prefix
        outputs = os.path.join(prefix, f"output_{name}/")
        inputs = os.path.join(prefix, f"input_{name}/")

        for path in [base, outputs, inputs]:
            os.makedirs(path, exist_ok=True)

        #os.path.
        base_prefix = os.path.join(base,name)
        outputs_prefix =  os.path.join(outputs,name)
        inputs_prefix =  os.path.join(inputs,name)

        return base_prefix, outputs_prefix, inputs_prefix


    def _decimate_name(self, input_name : str) -> str:
        name = input_name
        # Reduce the name size if its too verbose.
        if len(name) > 25:
            tokens = re.split(r'[_\-\s]', name)

            def sum_size(in_tokens):
                total_length = 0
                for t in in_tokens:
                    total_length += len(t)
                return total_length
            
            while(sum_size(tokens) > 25):
                reduced = sum_size(tokens)
                for i, t in enumerate(tokens):
                    if len(t) > 3:
                        tokens[i] = t[:-1]
                # Can't get smaller!
                if reduced == sum_size(tokens):
                    break

            tokens = [token for token in tokens if (token)]
            name = "_".join(tokens)

        return name

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
        return True# 
