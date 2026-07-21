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

from commec.config.constants import DEFAULT_CONFIG_YAML_PATH, DEPRECATED_CONFIG_KEYS
from commec.utils.file_utils import expand_and_normalize
from commec.utils.dict_utils import deep_update

logger = logging.getLogger(__name__)


class YamlIOValidationError(ValueError):
    """Custom exception for errors when handling input and output with Yaml Config."""


def load_config_from_yaml(config_filepath: str | os.PathLike) -> dict:
    """
    Loads a yaml file, ensuring it's a dictionary.
    """
    try:
        with open(config_filepath, "r", encoding="utf-8") as file:
            config_from_yaml = yaml.safe_load(file)
    except ParserError as e:
        raise ValueError(
            f"Invalid yaml syntax in configuration file: {config_filepath}"
        ) from e

    if not isinstance(config_from_yaml, dict):
        raise TypeError(
            f"Loaded configuration file did not result in a dictionary: {file}"
        )
    return config_from_yaml


def warn_and_strip_deprecated_keys(
    config_from_yaml: dict, config_filepath: str | os.PathLike
):
    """
    Remove any retired top-level keys from a loaded config, warning (rather than erroring)
    for each one found. This keeps older configs runnable across a deprecation window: the
    deprecated key is ignored instead of being rejected as an unrecognized key.
    """
    for key, reason in DEPRECATED_CONFIG_KEYS.items():
        if key in config_from_yaml:
            logger.warning(
                "Ignoring deprecated config key %r in %s: %s.",
                key,
                config_filepath,
                reason,
            )
            del config_from_yaml[key]


def update_config_from_yaml(
    existing_config: dict, config_filepath: str | os.PathLike
) -> dict:
    """
    Override YAML configuration based on provided YAML file. Items in the provided file, but
    not in the default YAML, will be ignored.
    """
    config_from_yaml = load_config_from_yaml(config_filepath)
    warn_and_strip_deprecated_keys(config_from_yaml, config_filepath)
    updated_config, rejected = deep_update(existing_config, config_from_yaml)
    if rejected:
        keys = ", ".join(f"{k}={v!r}" for k, v in rejected)
        raise YamlIOValidationError(
            f"Unrecognized key(s) in {config_filepath}: {keys}. "
            "Check for typos against the packaged default config."
        )
    return updated_config


def update_config_from_cli(yaml_config: dict, args: argparse.Namespace):
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
            logger.debug(
                "Command line arguments updated %s to: %s", arg, getattr(args, arg)
            )
            yaml_config[arg] = getattr(args, arg)
    return yaml_config


def format_config_paths(yaml_config: dict) -> dict:
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

    if not yaml_config.get("base_paths").get("default"):
        raise YamlIOValidationError(
            "base_paths.default is not configured. Set it via one of:\n"
            "  - `commec screen -d /path/to/databases ...`\n"
            "  - `base_paths: {default: /abs/path/to/databases}` in your --config YAML\n"
            "  - run `commec setup` to download databases and emit a config snippet"
        )
    try:

        def recursive_format(nested_yaml, base_paths):
            """
            Recursively apply string formatting to
            read paths from nested yaml config dicts.
            """
            if isinstance(nested_yaml, dict):
                return {
                    key: recursive_format(value, base_paths)
                    for key, value in nested_yaml.items()
                }
            if isinstance(nested_yaml, str):
                try:
                    return nested_yaml.format(**base_paths)
                except KeyError as e:
                    raise YamlIOValidationError(
                        f"Unknown base path key {e} referenced in path: {nested_yaml!r}. "
                        f"Known keys: {sorted(base_paths.keys())}"
                    ) from e
            return nested_yaml

        # Update the base paths
        yaml_config["base_paths"] = recursive_format(
            yaml_config["base_paths"], yaml_config["base_paths"]
        )

        # Every resolved base path must be absolute. Normalize with trailing separator.
        for key, value in yaml_config["base_paths"].items():
            if value is None:
                continue
            expanded = expand_and_normalize(value)
            if not os.path.isabs(expanded):
                raise YamlIOValidationError(
                    f"base_paths.{key} must be an absolute path, got: {value!r}"
                )
            yaml_config["base_paths"][key] = os.path.join(expanded, "")

        # Recursively format all paths
        yaml_config = recursive_format(yaml_config, yaml_config["base_paths"])

    except TypeError as e:
        logger.error(
            "Encountered unexpected TypeError during yaml config base path substitution: %s",
            e,
        )
        raise YamlIOValidationError(
            f"Encountered unexpected TypeError during yaml config base path substitution: {e}"
        )

    _validate_database_paths(yaml_config.get("databases", {}))

    return yaml_config


def _validate_database_paths(tree, breadcrumbs=()):
    """
    Walk the resolved `databases:` subtree and raise on any string value that
    isn't an absolute path.
    """
    if isinstance(tree, dict):
        for key, value in tree.items():
            _validate_database_paths(value, breadcrumbs + (key,))
        return
    if isinstance(tree, str) and breadcrumbs[-1] == "path":
        location = ".".join(breadcrumbs)
        if not os.path.isabs(expand_and_normalize(tree)):
            raise YamlIOValidationError(
                f"databases.{location} must be an absolute path, got: {tree!r}. "
                "Use absolute paths in YAML; CLI `-d` accepts relative paths."
            )


def get_defaults() -> dict:
    """
    Returns a deep copy of the defaults that can then be modified as needed.
    """
    default_yaml = importlib.resources.files("commec").joinpath(
        DEFAULT_CONFIG_YAML_PATH
    )
    default_config = None
    if default_yaml.exists():
        default_config = load_config_from_yaml(str(default_yaml))
        return default_config
    raise FileNotFoundError(
        f"No default yaml found. Expected at {DEFAULT_CONFIG_YAML_PATH}"
    )


if __name__ == "__main__":
    print("Commec yaml defaults:\n", pformat(get_defaults()))
