"""
Regression tests for the database download base URL, and for how `-d` and `-y`
divide the work of saying where the databases are: that an unusable URL is
rejected rather than crashing mid-download, and that a config is applied when it
is named and left alone when it is merely sitting in a directory that was.
"""

import argparse
import json
import os
from unittest.mock import patch

import pytest
import yaml

from commec.cli import ScreenArgumentParser
from commec.cli import main as commec_main
from commec.config.screen_io import ScreenIO
from commec.config.yaml_io import (
    YamlIOValidationError,
    default_base_url,
    format_config_paths,
    validate_base_url,
)
from commec.screen import add_args as screen_add_args
from commec.setup import (
    CommecDatabaseUpdater,
    CommecSetup,
    add_args,
    check_for_updates,
    resolve_base_url,
    unwritable_download_directories,
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


def test_setup_records_base_url_in_written_config(database_dir):
    """The config written by setup carries the URL the databases came from, at the
    top of the file, so a later run given it with -y finds the same host."""
    with patch("commec.setup.fetch_latest_revisions", return_value=LATEST) as fetched:
        CommecSetup(parse_setup_args(["-d", str(database_dir), "--base-url", MIRROR]))

    fetched.assert_called_once_with(MIRROR)
    written = yaml.safe_load((database_dir / "config.yaml").read_text())
    assert written["base_url"] == MIRROR
    assert next(iter(written)) == "base_url"


def test_setup_uses_the_url_from_the_config_it_is_given(database_dir, capsys):
    """A config, once named, is where the host comes from, which is how a mirror
    install keeps fetching from its mirror. A URL differing only by a trailing
    slash counts as recorded, and must not be nagged about."""
    config_path = database_dir / "config.yaml"
    config_path.write_text(
        f"base_url: {MIRROR}/\nbase_paths:\n  default: {database_dir}\n"
    )

    with patch("commec.setup.fetch_latest_revisions", return_value=LATEST) as fetched:
        setup = CommecSetup(parse_setup_args(["-y", str(config_path)]))

    assert setup.base_url == MIRROR
    fetched.assert_called_once_with(MIRROR)
    assert "does not record the base URL" not in capsys.readouterr().out


def test_setup_d_does_not_read_the_databases_directory_config(database_dir):
    """`-d` names a directory of databases and nothing more. A config left there
    may describe a different install, down to paths in another directory, so it is
    not read: fetching from the wrong host is recoverable, downloading over
    somebody else's databases is not."""
    (database_dir / "config.yaml").write_text(yaml.safe_dump({"base_url": MIRROR}))

    with patch("commec.setup.fetch_latest_revisions", return_value=LATEST) as fetched:
        setup = CommecSetup(parse_setup_args(["-d", str(database_dir)]))

    assert setup.base_url == DEFAULT
    fetched.assert_called_once_with(DEFAULT)


def test_setup_d_is_not_blocked_by_a_broken_directory_config(database_dir):
    """A stale, foreign, or unparseable config left in a databases directory must
    not stop a run that never asked for it. Refusing to update because of a file
    nobody pointed at is the failure this replaced."""
    (database_dir / "config.yaml").write_text("base_paths:\n\tdefault: /x\n")

    with patch("commec.setup.fetch_latest_revisions", return_value=LATEST):
        setup = CommecSetup(parse_setup_args(["-d", str(database_dir), "-m"]))

    assert setup.base_url == DEFAULT


def test_screen_d_does_not_read_the_databases_directory_config(database_dir):
    """`commec screen -d dir` screens against the databases in the directory it was
    given. A config found there may point at another directory entirely, so it is
    not read, and cannot quietly turn settings on either."""
    (database_dir / "config.yaml").write_text(
        yaml.safe_dump({"base_url": MIRROR, "auto_update_databases": True})
    )

    parser = ScreenArgumentParser()
    screen_add_args(parser)
    args = parser.parse_args([INPUT_QUERY, "-d", str(database_dir)])

    config = ScreenIO(args).config
    assert config["base_url"] == DEFAULT
    assert config["auto_update_databases"] is False


def test_d_and_y_combine_with_d_winning_the_base_path(database_dir):
    """The two flags are not mutually exclusive: the config says how the databases
    are laid out, -d says which directory this machine keeps them in. A path the
    config states outright is left alone, which is how a large database is kept on
    a separate volume -- the one thing -y can express that -d cannot."""
    elsewhere = database_dir / "bulk" / "nucleotide" / "nucl"
    config_path = database_dir / "split.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "base_paths": {"default": str(database_dir / "unused")},
                "databases": {"best_match": {"nucleotide": {"path": str(elsewhere)}}},
            }
        )
    )

    with patch("commec.setup.fetch_latest_revisions", return_value=LATEST):
        config = CommecSetup(
            parse_setup_args(["-d", str(database_dir), "-y", str(config_path), "-m"])
        ).config

    assert config["base_paths"]["default"] == os.path.join(str(database_dir), "")
    assert config["databases"]["best_match"]["nucleotide"]["path"] == str(elsewhere)
    assert config["databases"]["biorisk"]["path"] == str(
        database_dir / "biorisk" / "biorisk.hmm"
    )


def test_setup_requires_a_destination(capsys):
    """Neither flag leaves nowhere to download to, and no safe default to invent.
    This path was unreachable while the two flags were mutually exclusive, and its
    message was never formatted, so it is checked as rendered output."""
    with pytest.raises(SystemExit) as exit_info:
        CommecSetup(parse_setup_args([]))

    assert exit_info.value.code == 1
    printed = capsys.readouterr().out
    assert "Neither a databases directory" in printed
    assert "{ERROR_CHECK}" not in printed


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
    """Every way a -y config can be unusable is a user error, so each exits with
    one line rather than a traceback: a path that isn't there, one that cannot be
    read, and one that will not parse."""
    config_path = tmp_path / "config.yaml"
    if kind == "directory":
        config_path.mkdir()
    elif kind == "unparseable":
        config_path.write_text("base_paths:\n\tdefault: /x\n")

    with patch("sys.argv", ["commec", "setup", "-y", str(config_path), "-m"]):
        with pytest.raises(SystemExit) as exit_info:
            commec_main()

    assert exit_info.value.code == 1
    assert f"Configuration error: {message}" in capsys.readouterr().err


@pytest.mark.skipif(os.geteuid() == 0, reason="root ignores directory permissions")
def test_unwritable_download_directories_checks_each_database(tmp_path):
    """Checked per database, not just for the parent: a config may put one database
    on another volume, and a multi-GB download that only finds out at extraction
    time has already spent the transfer."""
    readonly = tmp_path / "readonly"
    readonly.mkdir()
    readonly.chmod(0o500)
    try:
        updaters = {
            "biorisk": CommecDatabaseUpdater(str(tmp_path / "dbs/biorisk"), DEFAULT),
            "best_match": CommecDatabaseUpdater(str(readonly / "best_match"), DEFAULT),
        }
        assert unwritable_download_directories(updaters) == [
            str(readonly / "best_match")
        ]
    finally:
        readonly.chmod(0o700)


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
