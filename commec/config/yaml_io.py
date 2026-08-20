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

import argparse
import importlib.resources
import logging
import os
from pprint import pformat
from urllib.parse import urlparse, urlunparse

import yaml

from commec.config.constants import (
    DEFAULT_CONFIG_YAML_PATH,
    DEPRECATED_CONFIG_KEYS,
    SETUP_EXAMPLE_CONFIG_FILENAME,
)
from commec.utils.dict_utils import deep_update
from commec.utils.file_utils import expand_and_normalize

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
    # Every parse failure, not just ParserError: a tab indent raises ScannerError,
    # an undefined alias ComposerError, and so on. None of those are ValueErrors,
    # so leaving them unconverted surfaces raw PyYAML tracebacks to users.
    except yaml.YAMLError as e:
        raise YamlIOValidationError(
            f"Invalid yaml syntax in configuration file {config_filepath}: {e}"
        ) from e
    # Unreadable rather than unparseable: no permission to read it, or a path that
    # named a directory. Still something the user can fix, so it is reported the
    # same one-line way as a syntax error rather than as a traceback.
    except OSError as e:
        raise YamlIOValidationError(
            f"Could not read configuration file {config_filepath}: {e}"
        ) from e

    # An empty file states nothing, so it is treated as stating nothing rather
    # than as broken. Chiefly so a setup killed part-way through writing a config
    # doesn't leave a zero-byte file that blocks the command used to repair it.
    if config_from_yaml is None:
        return {}

    if not isinstance(config_from_yaml, dict):
        raise YamlIOValidationError(
            f"Configuration file is not a set of settings: {config_filepath}"
        )
    return config_from_yaml


def note_unread_directory_config(
    database_dir: str | os.PathLike, config_filepath: str | os.PathLike | None = None
) -> None:
    """
    Note the example config `commec setup` leaves in a databases directory.

    `-d` names a directory of databases and nothing more. The file setup writes
    there is an example describing that install, rewritten by every setup run, so
    it is not read: settings are applied by copying it and naming the copy with `-y`.

    Logged rather than printed, so a run that never wanted the file stays quiet
    while the reason is still in the log for a run that expected it to apply.

    Silent when that file is the one `-y` named, which is not unread at all.
    """
    directory_config = os.path.join(database_dir, SETUP_EXAMPLE_CONFIG_FILENAME)
    if not os.path.isfile(directory_config):
        return
    if config_filepath and os.path.abspath(config_filepath) == os.path.abspath(
        directory_config
    ):
        return
    logger.info(
        "Not reading %s; it is an example rewritten by every `commec setup` run."
        " Copy it, edit the copy, and pass that with -y/--config to apply its settings",
        directory_config,
    )


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
    # Asking for a config that turns out to be empty is always a mistake, and
    # saying so beats the puzzling "base_paths.default is not configured" that
    # would otherwise follow. Callers that discover a config rather than being
    # pointed at one check for emptiness themselves and never get here.
    if not config_from_yaml:
        raise YamlIOValidationError(
            f"{config_filepath} is empty: it contains no settings."
        )
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


def validate_base_url(base_url, source: str = "base_url") -> str:
    """
    Check that a database download base URL is one urllib can actually open, and
    return it as urllib will parse it, without a trailing slash.

    A URL urllib rejects does not fail politely: a missing scheme raises a bare
    ValueError and a bad port raises http.client.InvalidURL, neither of which is
    a URLError or an OSError, so neither is caught by any of the download error
    handling in `commec setup`. Anything unusable has to be rejected up front.

    `source` names where the value came from, so the message can point at the
    command line flag, the config key, or a specific file.
    """
    if not isinstance(base_url, str):
        raise YamlIOValidationError(
            f"{source} must be a URL string, got {type(base_url).__name__}: {base_url!r}"
        )
    if not base_url.strip():
        raise YamlIOValidationError(
            f"{source} is empty. Give a URL such as https://databases.commec.io, "
            "or leave it unset to use the default."
        )

    stripped = base_url.strip()

    # urlsplit deletes tabs and newlines from anywhere in the string, so a URL
    # holding one parses as a different URL than the one written down. Reject it
    # rather than silently using a host nobody typed.
    if any(character in stripped for character in "\t\r\n"):
        raise YamlIOValidationError(
            f"{source} contains a tab or newline, got: {base_url!r}"
        )

    parsed = urlparse(stripped)
    if parsed.scheme not in ("http", "https"):
        raise YamlIOValidationError(
            f"{source} must start with http:// or https://, got: {base_url!r}"
        )
    if not parsed.hostname:
        raise YamlIOValidationError(f"{source} has no host, got: {base_url!r}")
    try:
        parsed.port  # Accessing it is the check: a non-numeric port raises here.
    except ValueError as e:
        raise YamlIOValidationError(
            f"{source} has an invalid port, got: {base_url!r}"
        ) from e
    if parsed.username or parsed.password:
        raise YamlIOValidationError(
            f"{source} must not embed credentials: urllib leaves them in the request"
            f" host so the download fails, and they would be stored in plain text."
            f" Got: {base_url!r}"
        )
    if parsed.query or parsed.fragment or parsed.params:
        raise YamlIOValidationError(
            f"{source} must be a plain host and path, with no query or fragment,"
            f" got: {base_url!r}"
        )

    # Return what urllib will parse, rather than the raw input, so the value that
    # is compared, stored and printed is the value that gets used.
    return urlunparse(parsed).rstrip("/")


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

    An optional `base_url` is deliberately left out of the substitution: it is a
    URL rather than a path, and is validated by whatever is about to download
    from it, so a bad one cannot block commands that never download.
    """

    if not yaml_config.get("base_paths").get("default"):
        raise YamlIOValidationError(
            "base_paths.default is not configured. Set it via one of:\n"
            "  - `commec screen -d /path/to/databases ...`\n"
            "  - `base_paths: {default: /abs/path/to/databases}` in your --config YAML\n"
            "  - run `commec setup` to download databases and emit a config snippet"
        )

    try:

        def recursive_format(nested_yaml, base_paths, skip=()):
            """
            Recursively apply string formatting to
            read paths from nested yaml config dicts.

            Keys in `skip` are passed through untouched, for values that are not
            paths and must not have base paths substituted into them.
            """
            if isinstance(nested_yaml, dict):
                return {
                    key: value if key in skip else recursive_format(value, base_paths)
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

        # Recursively format all paths. base_url is skipped: it is a URL rather
        # than a path, so substituting base paths into it would splice a local
        # directory into a download address, and any brace it contains would
        # raise a formatting error naming neither the key nor the file.
        yaml_config = recursive_format(
            yaml_config, yaml_config["base_paths"], skip=("base_url",)
        )

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


def default_base_url() -> str:
    """
    The host databases are fetched from when nothing overrides it.

    Read from the packaged config rather than held as a constant, so that file
    is the single visible place the default is written down. For configs that
    were built by hand rather than merged from the defaults, and so carry no
    base_url of their own.
    """
    return get_defaults()["base_url"]


if __name__ == "__main__":
    print("Commec yaml defaults:\n", pformat(get_defaults()))
