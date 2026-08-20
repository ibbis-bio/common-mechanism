"""
Regression tests for the database download base URL, and for the example config
`commec setup` writes into the databases directory: that an unusable URL is
rejected rather than crashing mid-download, that the example config is rewritten
by every run and records the host actually used, and that a config sitting in a
databases directory is never read by the commands that were merely pointed at
that directory.
"""

import argparse
import json
import os
from unittest.mock import patch

import pytest
import yaml

from commec.cli import ScreenArgumentParser
from commec.cli import main as commec_main
from commec.config.constants import SETUP_EXAMPLE_CONFIG_FILENAME
from commec.config.screen_io import ScreenIO
from commec.config.yaml_io import (
    YamlIOValidationError,
    default_base_url,
    format_config_paths,
    validate_base_url,
)
from commec.screen import add_args as screen_add_args
from commec.setup import (
    CommecSetup,
    add_args,
    check_for_updates,
    resolve_base_url,
)

DEFAULT = default_base_url()
MIRROR = "https://mirror.example.org/commec"
LATEST = {"latest": {"biorisk": "1.0"}}
INPUT_QUERY = os.path.join(os.path.dirname(__file__), "test_data/single_record.fasta")


def parse_setup_args(argv):
    parser = argparse.ArgumentParser()
    add_args(parser)
    return parser.parse_args(argv)


@pytest.fixture
def database_dir(tmp_path):
    """A database directory holding one up-to-date database, so a setup run
    completes without downloading anything."""
    biorisk = tmp_path / "biorisk"
    biorisk.mkdir()
    (biorisk / "manifest.json").write_text(
        json.dumps({"component": "biorisk", "revision": "1.0"})
    )
    return tmp_path


@pytest.mark.parametrize(
    "value, message",
    [
        # urlopen raises a bare ValueError for these, which no download error
        # handling catches, so they have to be rejected up front.
        ("mirror.example.org", "must start with http:// or https://"),
        ("https://", "has no host"),
        # An empty value is usually an unset variable, not a request for the default.
        ("", "is empty"),
        # urlopen raises http.client.InvalidURL for these, which is neither an
        # OSError nor a URLError, so the download paths cannot catch it either.
        ("https://host:notaport", "invalid port"),
        ("https://user:pass@host", "must not embed credentials"),
        # urlsplit deletes these mid-string, so the URL parsed is not the one written.
        ("https://evil.example\n.commec.io", "contains a tab or newline"),
    ],
)
def test_validate_base_url_rejects(value, message):
    with pytest.raises(YamlIOValidationError, match=message):
        validate_base_url(value)


def test_resolve_base_url_precedence():
    """Command line beats the config, which beats the built-in default. Only an
    absent value falls through; anything present is used or rejected."""
    assert resolve_base_url(MIRROR, {"base_url": "https://from-config"}) == MIRROR
    assert resolve_base_url(None, {"base_url": MIRROR}) == MIRROR
    assert resolve_base_url(None, {"base_url": None}) == DEFAULT
    assert resolve_base_url(None, {}) == DEFAULT


def test_base_url_is_not_treated_as_a_path():
    """A URL must survive the base_paths substitution pass untouched, braces and
    all, rather than having a local directory spliced into a download address. It
    is also not validated here, so a bad one cannot block `screen` or `list`,
    which never download anything."""
    for url in ("mirror.local", "https://m.example.org/{default}/x"):
        config = {"base_paths": {"default": "/commec-dbs/"}, "base_url": url}
        assert format_config_paths(config)["base_url"] == url


def test_setup_writes_an_example_config_and_says_to_copy_it(database_dir, capsys):
    """The config setup leaves behind describes this install: the databases
    directory and the URL they came from, at the top of the file. Nothing in the
    file says it is a copy to work from -- the run says so instead -- so the advice
    is checked as printed output."""
    with patch("commec.setup.fetch_latest_revisions", return_value=LATEST) as fetched:
        CommecSetup(parse_setup_args(["-d", str(database_dir), "--base-url", MIRROR]))

    fetched.assert_called_once_with(MIRROR)
    example_config = database_dir / SETUP_EXAMPLE_CONFIG_FILENAME
    assert example_config.name == "example-config.yaml"

    printed = capsys.readouterr().out
    assert "Every setup run rewrites it" in printed
    assert f"cp {example_config} my-config.yaml" in printed

    written = yaml.safe_load(example_config.read_text())
    assert written["base_url"] == MIRROR
    assert next(iter(written)) == "base_url"
    assert written["base_paths"]["default"] == str(database_dir.resolve())


def test_setup_rewrites_the_example_config_every_run(database_dir):
    """Rewritten on an update run, and over a file already there, so it never
    describes an install other than the one in this directory. The previous
    behaviour -- refusing to overwrite -- left a stale config in place for as long
    as the directory survived."""
    example_config = database_dir / SETUP_EXAMPLE_CONFIG_FILENAME
    example_config.write_text("base_url: https://stale.example.org\n")

    with patch("commec.setup.fetch_latest_revisions", return_value=LATEST):
        # Nothing to download: the fixture's database is already up to date, so
        # this is the update-check run that used to write nothing at all.
        CommecSetup(parse_setup_args(["-d", str(database_dir), "--base-url", MIRROR]))

    assert yaml.safe_load(example_config.read_text())["base_url"] == MIRROR

    with patch("commec.setup.fetch_latest_revisions", return_value=LATEST):
        CommecSetup(parse_setup_args(["-d", str(database_dir)]))

    assert yaml.safe_load(example_config.read_text())["base_url"] == DEFAULT


def test_setup_mock_run_writes_nothing(database_dir):
    """`-m/--mock` changes nothing on disk, the example config included, so it can
    be used to preview an update against a directory in use."""
    with patch("commec.setup.fetch_latest_revisions", return_value=LATEST):
        CommecSetup(parse_setup_args(["-d", str(database_dir), "-m"]))

    assert not (database_dir / SETUP_EXAMPLE_CONFIG_FILENAME).exists()


def test_setup_does_not_read_the_databases_directory_config(database_dir):
    """`-d` names a directory of databases and nothing more. The config left there
    is an example for the user to copy, not an input: reading it back would let a
    stale or hand-edited file redirect the next download."""
    (database_dir / SETUP_EXAMPLE_CONFIG_FILENAME).write_text(
        yaml.safe_dump({"base_url": MIRROR})
    )

    with patch("commec.setup.fetch_latest_revisions", return_value=LATEST) as fetched:
        setup = CommecSetup(parse_setup_args(["-d", str(database_dir)]))

    assert setup.base_url == DEFAULT
    fetched.assert_called_once_with(DEFAULT)


def test_setup_is_not_blocked_by_a_broken_directory_config(database_dir):
    """A stale, foreign, or unparseable config left in a databases directory must
    not stop a run that never asked for it. Refusing to update because of a file
    nobody pointed at is the failure this replaced."""
    (database_dir / SETUP_EXAMPLE_CONFIG_FILENAME).write_text(
        "base_paths:\n\tdefault: /x\n"
    )

    with patch("commec.setup.fetch_latest_revisions", return_value=LATEST):
        setup = CommecSetup(parse_setup_args(["-d", str(database_dir), "-m"]))

    assert setup.base_url == DEFAULT


def test_screen_does_not_read_the_databases_directory_config(database_dir):
    """`commec screen -d dir` screens against the databases in the directory it was
    given. The example config found there is a reference copy, not settings: it
    cannot quietly turn features on or point at another host."""
    (database_dir / SETUP_EXAMPLE_CONFIG_FILENAME).write_text(
        yaml.safe_dump({"base_url": MIRROR, "auto_update_databases": True})
    )

    parser = ScreenArgumentParser()
    screen_add_args(parser)
    args = parser.parse_args([INPUT_QUERY, "-d", str(database_dir)])

    config = ScreenIO(args).config
    assert config["base_url"] == DEFAULT
    assert config["auto_update_databases"] is False


def test_setup_requires_a_databases_directory(capsys):
    """Setup downloads a fixed set of databases into one parent directory, and has
    no safe default to invent for where that is, so -d is not optional."""
    with pytest.raises(SystemExit) as exit_info:
        parse_setup_args([])

    assert exit_info.value.code == 2
    assert "-d/--databases" in capsys.readouterr().err


def test_setup_does_not_take_a_config(database_dir, capsys):
    """Setup no longer reads a config at all: it installs the packaged layout under
    -d. A -y that was silently ignored would be worse than one that is rejected."""
    with pytest.raises(SystemExit) as exit_info:
        parse_setup_args(["-d", str(database_dir), "-y", "config.yaml"])

    assert exit_info.value.code == 2
    assert "unrecognized arguments: -y" in capsys.readouterr().err


@pytest.mark.parametrize(
    "kind, message",
    [
        ("missing", "--config YAML not found"),
        ("directory", "Could not read"),
        ("unparseable", "Invalid yaml syntax"),
    ],
)
def test_cli_reports_an_unusable_config_without_a_traceback(
    tmp_path, capsys, kind, message
):
    """Every way a -y config given to `screen` can be unusable is a user error, so
    each exits with one line rather than a traceback: a path that isn't there, one
    that cannot be read, and one that will not parse."""
    config_path = tmp_path / "my-config.yaml"
    if kind == "directory":
        config_path.mkdir()
    elif kind == "unparseable":
        config_path.write_text("base_paths:\n\tdefault: /x\n")

    argv = [
        "commec",
        "screen",
        "-y",
        str(config_path),
        "-o",
        str(tmp_path / "out"),
        INPUT_QUERY,
    ]
    with patch("sys.argv", argv):
        with pytest.raises(SystemExit) as exit_info:
            commec_main()

    assert exit_info.value.code == 1
    assert f"Configuration error: {message}" in capsys.readouterr().err


def test_check_for_updates_uses_configured_base_url(database_dir):
    """The path `commec screen` takes for automatic updates must honour the
    recorded mirror, not go back to the public default."""
    config = {
        "base_url": MIRROR,
        "databases": {"biorisk": {"path": str(database_dir / "biorisk/biorisk.hmm")}},
    }

    with patch("commec.setup.fetch_latest_revisions", return_value=LATEST) as fetched:
        _, updaters = check_for_updates(config)

    fetched.assert_called_once_with(MIRROR)
    assert updaters["biorisk"].base_url == MIRROR
