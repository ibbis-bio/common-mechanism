#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Defines helper functions for the ingestion of the config yaml file, needed
for commec screen as well as commec list.
The yaml config file import involves the following features:
- default values as determined by the screen-default-config.yaml
- Ability to selectively override defaults
- Ability to parse additional path substitution
"""
import os
from pprint import pformat
import logging
import argparse
import importlib.resources
import yaml
from yaml.parser import ParserError

from commec.config.constants import DEFAULT_CONFIG_YAML_PATH
from commec.utils.dict_utils import deep_update

logger = logging.getLogger(__name__)

def load_config_from_yaml(config_filepath: str | os.PathLike) -> dict:
    """
    Loads a yaml file, ensuring it's a dictionary.
    """
    try:
        with open(config_filepath, "r", encoding = "utf-8") as file:
            config_from_yaml = yaml.safe_load(file)
    except ParserError as e:
        raise ValueError(f"Invalid yaml syntax in configuration file: {config_filepath}") from e

    if not isinstance(config_from_yaml, dict):
        raise TypeError(f"Loaded configuration file did not result in a dictionary: {file}")
    return config_from_yaml

def update_config_from_yaml(existing_config : dict, config_filepath: str | os.PathLike) -> dict:
    """
    Override YAML configuration based on provided YAML file. Items in the provided file, but
    not in the default YAML, will be ignored.
    """
    config_from_yaml = load_config_from_yaml(config_filepath)
    updated_config, rejected = deep_update(existing_config, config_from_yaml)
    for rejects in rejected:
        logger.warning("The follow input from the user provided"
            " configuration was not recognised: %s : %s",
            rejects[0], rejects[1])
    return updated_config

def update_config_from_cli(yaml_config : dict, args: argparse.Namespace):
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
        if arg in yaml_config and hasattr(args, arg):
            logger.debug("Command line arguments updated %s to: %s", arg, getattr(args,arg))
            yaml_config[arg] = getattr(args, arg)
    return yaml_config

def format_config_paths(yaml_config : dict) -> dict:
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

    if not yaml_config.get("base_paths"):
        logger.debug("No Base paths to perform yaml substitution.")
        return yaml_config
    
    try:
        base_paths = yaml_config["base_paths"]

        # Ensure all the base paths end with a separator "/""
        for key, value in base_paths.items():
            base_paths[key] = os.path.join(value,'')

        def recursive_format(nested_yaml, base_paths):
            """
            Recursively apply string formatting to
            read paths from nested yaml config dicts.
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

        # Recursively format all paths
        yaml_config["base_paths"] = recursive_format(yaml_config["base_paths"], yaml_config["base_paths"])
        yaml_config = recursive_format(yaml_config, yaml_config["base_paths"])

        yaml_config = recursive_format(yaml_config, base_paths)
    except TypeError as e:
        logger.error("Encountered unexpected TypeError during yaml config base path substitution: %s", e)
        pass

    return yaml_config

def get_defaults() -> dict:
    """
    Returns a deep copy of the defaults that can then be modified as needed.
    """
    default_yaml = importlib.resources.files("commec").joinpath(DEFAULT_CONFIG_YAML_PATH)
    default_config = None
    if default_yaml.exists():
        default_config = load_config_from_yaml(str(default_yaml))
        return default_config
    raise FileNotFoundError(
        f"No default yaml found. Expected at {DEFAULT_CONFIG_YAML_PATH}"
        )

if __name__ == "__main__":
    print("Commec yaml defaults:\n", pformat(get_defaults()))
