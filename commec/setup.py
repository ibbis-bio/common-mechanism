#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Module for CLI setup of Commec, such that required
databases are downloaded in a desired database directory.
"""

import argparse
import hashlib
import importlib.resources
import json
import os
import shutil
import ssl
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from urllib import error, request
from urllib.parse import urlparse

import yaml

from commec import __version__ as VERSION
from commec.config import yaml_io as YamlIO
from commec.config.constants import DEFAULT_CONFIG_YAML_PATH, SETUP_CONFIG_FILENAME
from commec.config.yaml_io import YamlIOValidationError
from commec.utils.file_utils import expand_and_normalize, remove_filename_from_path

DESCRIPTION = """Helper script for downloading or updating the databases
 required for running the Common Mechanism Screen"""

C_F_ORANGE = "\033[38;5;202m"  # Colour Foreground Orange.
C_F_AMBER = "\033[38;5;214m"  # Colour Foreground Amber
C_F_GRAY = "\033[38;5;242m"  # Colour Foreground Gray
C_F_BLUE = "\033[38;5;17m"  # Colour Foreground Blue
C_B_BLUE = "\033[48;5;17m"  # Colour Background Blue
C_RESET = "\033[0m"  # Reset Console Formatting.
C_BOLD = "\033[1m"  # Reset Console Formatting.

# Something failed: the run is over, or this database won't be updated.
ERROR_CHECK = C_F_ORANGE + C_BOLD + " X " + C_RESET
# Something is worth knowing but the run carries on regardless. Kept visually
# distinct so a successful run needing a nudge doesn't read as a failed one.
WARNING_CHECK = C_F_AMBER + C_BOLD + " ! " + C_RESET
STEP = " ➔ "
BULLET = "  ● "

SPLASH_IMAGE = f"""
                        Welcome to

     ██████╗ ██████╗ ███╗   ███╗███╗   ███╗███████╗ ██████╗ {C_F_ORANGE}         ▄▄               {C_RESET}
    ██╔════╝██╔═══██╗████╗ ████║████╗ ████║██╔════╝██╔════╝ {C_F_ORANGE}       ▄███▌              {C_RESET}
    ██║     ██║   ██║██╔████╔██║██╔████╔██║█████╗  ██║      {C_F_ORANGE}      ▐█████              {C_RESET}
    ██║     ██║   ██║██║╚██╔╝██║██║╚██╔╝██║██╔══╝  ██║      {C_F_ORANGE}     ▐██████▌             {C_RESET}
    ╚██████╗╚██████╔╝██║ ╚═╝ ██║██║ ╚═╝ ██║███████╗╚██████╗ {C_F_ORANGE}     ███████▌             {C_RESET}
     ╚═════╝ ╚═════╝ ╚═╝     ╚═╝╚═╝     ╚═╝╚══════╝ ╚═════╝ {C_F_ORANGE}    ▐███████▌             {C_RESET}
{C_B_BLUE}    █ █████▄ █████▄ █ ▄█▀█▄                                 {C_F_ORANGE}     ███████   ▄█▄      {C_RESET}
{C_B_BLUE}    █ █    █ █    █ █ █   ▀          DATABASE               {C_F_ORANGE}      █████▌  ▄███▄▄     {C_RESET}  
{C_B_BLUE}    █ █████▄ █████▄ █ ▀███▄            UPDATE               {C_F_ORANGE}      ▐█████▄██▀    ▀▄    {C_RESET}   
{C_B_BLUE}    █ █    █ █    █ █ ▄   █              UTILITY            {C_F_ORANGE}      ▐████████       ▌  {C_RESET}
{C_B_BLUE}    █ █████▀ █████▀ █ ▀█▄█▀                                 {C_F_ORANGE}      ████████▀         {C_RESET}
                                                            {C_F_ORANGE}   ▄▄██████▀▀             {C_RESET}
                                                            {C_F_ORANGE} ▀▀                       {C_RESET}"""

SPLASH_TEXT = """                    The Common Mechanism!
\nInternational Biosecurity and Biosafety Initiative for Science Copyright © 2021-2024
\nThis script will help download or update the mandatory databases required for using 
Commec Screen, and requires a stable internet connection.\n"""


def add_args(input_options: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """
    Add module arguments to an ArgumentParser object.
    """

    # Deliberately not a mutually exclusive group: a config can describe the
    # database layout while -d says which directory this machine keeps it in, so
    # one config held in version control serves hosts that put their databases in
    # different places. At least one of the two is still required, which argparse
    # cannot express for "either or both", so CommecSetup checks it instead.
    input_options.add_argument(
        "-d",
        "--databases",
        dest="database_dir",
        default=None,
        help="Path to a parent directory to download, or update, the Commec databases."
        " Overrides base_paths.default in any --config YAML.",
    )
    input_options.add_argument(
        "-y",
        "--config",
        dest="config_yaml",
        default="",
        help="Commec yaml configuration file describing where the databases live."
        " Combine with -d to keep the parent directory out of the file. Cannot"
        " change which databases are downloaded.",
    )

    # Optional
    input_options.add_argument(
        "-m",
        "--mock",
        dest="dry_run",
        default=False,
        action="store_true",
        help="Mock output potential update information, without performing any download or update.",
    )

    input_options.add_argument(
        "-e",
        "--experimental",
        dest="experimental",
        default=False,
        action="store_true",
        help="Download the latest experimental databases. Note: These are not intended for general use.",
    )

    input_options.add_argument(
        "-r",
        "--restore",
        dest="restore_json",
        default=None,
        help="Point setup to a commec screen output json file to replicate the databases used for that run.",
    )

    input_options.add_argument(
        "--base-url",
        dest="base_url",
        default=None,
        type=_base_url_arg,
        help="Base URL to download the databases.",
    )

    return input_options


def _base_url_arg(value: str) -> str:
    """
    argparse `type` for --base-url, so an unusable URL is rejected at parse time
    with the reason shown, rather than surfacing much later as a download failure.
    """
    try:
        return YamlIO.validate_base_url(value, "--base-url")
    except ValueError as e:
        raise argparse.ArgumentTypeError(str(e)) from e


class CommecSetup:
    def __init__(self, args):
        working_directory = Path()
        # Updaters dictionary, of {[name : UpdaterObject]}
        updaters = {}
        self.config = YamlIO.get_defaults()
        # The config named by -y, or None. A config sitting in the databases
        # directory is not read: -d names a directory of databases, and a file
        # found in one may describe another install, down to paths that would
        # send this download somewhere else entirely.
        self.config_filepath = None

        # Neither flag leaves no destination to download to, and there is no safe
        # default to invent for one.
        if not args.database_dir and not args.config_yaml:
            print(
                f"{ERROR_CHECK}Neither a databases directory (-d/--databases) nor a"
                " configuration yaml (-y/--config) was provided. Exiting now."
            )
            sys.exit(1)

        if args.config_yaml:
            self.config_filepath = args.config_yaml
            # -d wins for the parent directory, so a config may leave
            # base_paths.default unset and still be used on a machine that keeps
            # its databases elsewhere. Paths the config states outright are left
            # alone, which is how a database is kept on a separate volume.
            self.config = read_config_yaml_for_database_setup(
                args.config_yaml,
                database_dir=(
                    expand_and_normalize(args.database_dir)
                    if args.database_dir
                    else None
                ),
            )

        if args.database_dir:
            working_directory = expand_and_normalize(args.database_dir)
            YamlIO.note_unread_directory_config(working_directory, self.config_filepath)
        else:
            working_directory = expand_and_normalize(
                self.config.get("base_paths").get("default")
            )

        cli_base_url = getattr(args, "base_url", None)

        # The config that records this install, and so has to carry the base URL
        # for a later run given it with -y: the one named on the command line, or
        # the one a first-time setup will write into the databases directory.
        self.recording_config_path = self.config_filepath or os.path.join(
            working_directory, SETUP_CONFIG_FILENAME
        )
        # Whether that file was there before this run started. Only for telling
        # "was already there" from "appeared while this run was downloading" when
        # explaining a skipped write, since a download can take hours.
        self.config_present_at_start = os.path.exists(self.recording_config_path)
        self.recorded_base_url = normalized_base_url(
            self.config.get("base_url"), f"base_url in {self.recording_config_path}"
        )

        # Resolved after the config is read, so a base_url recorded in it is
        # picked up, and recorded back so any config written out carries it.
        try:
            self.base_url = resolve_base_url(
                cli_base_url, self.config, f"base_url in {self.recording_config_path}"
            )
        except ValueError as e:
            print(f"{ERROR_CHECK}{e}")
            sys.exit(1)
        self.config["base_url"] = self.base_url

        # Always, rather than only for a custom host: naming the host on every
        # run is what makes an install that has quietly gone back to the default
        # visible, and it costs one line.
        print(f"{STEP}Fetching databases from {self.base_url}")

        dry_run = getattr(args, "dry_run", False)

        if self.config_filepath:
            updaters = create_updaters_from_config(self.config, self.base_url)

        # Ensure working_directory is a valid path (even if it doesn't exist yet)
        if not check_directory_is_writable(working_directory):
            print(
                f"{ERROR_CHECK}Failed to parse working directory: {working_directory}"
            )

        databases = None
        # Fetch latest.json from R2.abs
        if args.restore_json:
            if not os.path.isfile(args.restore_json):
                print(
                    f"{ERROR_CHECK}Provided argument for import revisions from JSON file not a valid file: {args.restore_json}"
                )
                raise FileNotFoundError(args.restore_json)
            print(
                f"{STEP}A json database file was provided, fetching revision information..."
            )
            databases = fetch_revisions_from_json(args.restore_json)
            if not databases:
                print(
                    f"{ERROR_CHECK}Input JSON had no database revision information. Ensure the expected format: e.g.\n"
                    "{\n"
                    '  "database_info" : {\n'
                    '    "revisions" : {\n'
                    '       "biorisk" : "1.0",\n'
                    '       "best_match" : "1.0",\n'
                    '       "low_concern" : "1.0",\n'
                    '       "control_lists" : "1.0"\n'
                    "    }\n  }\n}"
                )
                return
        else:
            print(f"{STEP}Fetching latest commec database revision information ...")
            latest = fetch_latest_revisions(self.base_url)

            if not latest:
                print(
                    f"{ERROR_CHECK} Could not fetch latest revisions. Check internet connection and try again."
                )
                return

            # For each database, create an updater.
            if args.experimental:  # == "latest" or args.revision == "experimental":
                databases = latest["experimental"]
            else:
                databases = latest["latest"]

        if not databases:
            print(
                f"{ERROR_CHECK}An issue occured when fetching databases. Exiting now."
            )
            return

        for database_name, revision in databases.items():
            updater = updaters.get(database_name)
            if not updater:
                download_location = os.path.join(working_directory, database_name)
                updater = CommecDatabaseUpdater(download_location, self.base_url)
                updaters[database_name] = updater
                if self.config_filepath:
                    print(
                        f"{WARNING_CHECK}{database_name} not present in"
                        f" {self.config_filepath}, it may be out of date."
                    )
            updater.name = database_name
            updater.check_for_update(revision)

        # Check for orphaned databases:
        for db_name, updater in updaters.items():
            if db_name not in databases.keys():
                print(
                    f"{WARNING_CHECK}{db_name} has no formal update route:"
                    f" {self.config_filepath} has deprecated database entries."
                )
                updater.update_required = False

        # Every directory that will actually be written to, not just the parent:
        # a config may place one database on another volume, and discovering at
        # extraction time that it cannot be written has already cost the download.
        for directory in unwritable_download_directories(updaters):
            print(f"{ERROR_CHECK}Cannot write to database directory: {directory}")

        updated_required = False
        # Each updater, reports on its status, does it need to update?
        print(f"{STEP}The following actions have been identified ...")
        for _db_name, updater in updaters.items():
            print(BULLET, updater.update_message)
            updated_required = updated_required or updater.update_required

        # Log output of intended changes.
        if dry_run:
            self.advise_on_base_url()
            print(
                f"\n{C_F_GRAY}This was a mock run. Run {C_F_ORANGE}'commec"
                f" setup'{C_F_GRAY} again without {C_F_ORANGE}'-m/--mock'{C_F_GRAY}"
                f" to continue with the update.{C_RESET}"
            )
            return

        # Perform the database updates, consider making this async.
        os.makedirs(working_directory, exist_ok=True)
        for database_name, updater in updaters.items():
            updater.perform_update()

        # A setup that brought its own config with -y already has one recording
        # this install, so nothing is written; otherwise write one, recording what
        # this run resolved so it can be handed back with -y.
        if self.config_filepath is None:
            self.config["base_paths"]["default"] = str(
                Path(working_directory).resolve()
            )
            # Any file at all stops the write, empty or not: a config we did not
            # create is never written over. Checked again here rather than
            # trusting the look taken at the start of the run, since the
            # downloads in between can take hours and another setup into the
            # same directory may have created one since.
            if os.path.exists(self.recording_config_path):
                reason = (
                    "was already there"
                    if self.config_present_at_start
                    else "appeared while this run was downloading"
                )
                print(
                    f"{WARNING_CHECK}{self.recording_config_path} {reason}, so no"
                    " config was written."
                    "\n   Remove it and run setup again to have one created."
                )
            else:
                with open(
                    self.recording_config_path, "w", encoding="utf-8"
                ) as output_config_file:
                    # Unsorted, so the written config keeps the same key order as
                    # the packaged default, base_url included.
                    yaml.safe_dump(self.config, output_config_file, sort_keys=False)
                # Nothing reads this file on its own, so saying it exists is only
                # half the message: the invocation that uses it is the other half.
                print(
                    f"{STEP}A config.yaml recording this install was created at"
                    f" {self.recording_config_path}."
                    "\n   Nothing reads it automatically. To screen with these"
                    " settings, and to"
                    "\n   keep updates coming from the same host, pass it back:"
                    f"\n     commec screen -y {self.recording_config_path} <input.fasta>"
                )
                # Only once a config has actually been written: it carries the
                # resolved URL, so there is nothing left to advise about. When
                # the write was skipped the advice still needs to be given.
                self.recorded_base_url = self.base_url

        self.advise_on_base_url()

        print(f"{STEP}Update check complete! Have a Biosafe-and-secure day!")

    def advise_on_base_url(self) -> None:
        """
        Point at the config that needs to record the base URL, when this run
        isn't the one writing that config.

        Shared by real and mock runs, so `-m/--mock` previews the same advice
        rather than leaving it to be discovered on the run that matters.
        """
        # No config yet means setup is about to write one, or would have on a
        # real run, and a config setup writes already carries the URL.
        if not os.path.isfile(self.recording_config_path):
            return

        warn_base_url_not_recorded(
            self.recording_config_path, self.base_url, self.recorded_base_url
        )


def fetch_supported_revisions() -> dict | None:
    """
    The commec package has the database_compatibility.yaml, which
    holds a datastructure that can be interrogated for latest revision
    information supported by this version of commec.
    Revisions contains a major, and minor revision. By definition, the
    major revision increments only when a commec update is required.

    Returns a dictionary where each database key contains the tuple:
    (min, max) where max is non-inclusive. Fulfilling the ">=min, <max" logic.
    """
    data_filename = importlib.resources.files("commec").joinpath(
        "database_compatibility.yaml"
    )
    supported_revision_data = YamlIO.load_config_from_yaml(data_filename)
    for db_name, revision_string in supported_revision_data[
        "supported_database_revisions"
    ].items():
        # Expects the following str format: ">=1.0,<2.0"
        major_min = revision_string.split(",")[0].split(".")[0][2:] + ".0"
        major_max = revision_string.split(",")[1].split(".")[0][1:] + ".0"
        supported_revision_data[db_name] = (
            DatabaseRevision(major_min),
            DatabaseRevision(major_max),
        )
    return supported_revision_data


def fetch_revisions_from_json(filename: str | os.PathLike) -> dict | None:
    """
    Any json output from commec screen should contain the revision information of the databases used.
    Expected json format is:
    database_info : {
        revisions : {
            biorisk : "1.0",
            best_match : "1.0",
            low_concern : "1.0",
            control_list : "1.0"
        }
    }
    Returns an empty dict on failure.
    """
    json_string: str
    with open(filename, "r", encoding="utf-8") as json_file:
        json_string = json_file.read()
    my_data: dict = json.loads(json_string)
    db_revisions = my_data.get("database_info", {})
    db_revisions = db_revisions.get("revisions", {})
    return db_revisions


def fetch_latest_revisions(base_url: str) -> dict | None:
    """
    A latest.json manifest exists at the root of the R2 bucket, listing
    the latest revision of each database. Downloads and parses it into a
    dict. Returns None if it could not be fetched or parsed.
    """
    raw = fetch_r2_object("latest.json", base_url=base_url)
    if raw is None:
        return None

    try:
        return json.loads(raw.decode("utf-8"))
    except (json.JSONDecodeError, UnicodeDecodeError) as e:
        print(f"Invalid JSON in latest.json manifest: {e}")
        return None


def read_manifest(check_location):
    """
    All databases should have a manifest.json file.
    """
    manifest_location = remove_filename_from_path(check_location)
    manifest_filename = os.path.join(manifest_location, "manifest.json")
    manifest = {"component": "invalid", "revision": "0.0"}
    if os.path.exists(manifest_location) and os.path.isfile(manifest_filename):
        with open(manifest_filename, "r", encoding="utf-8") as manifest_file:
            json_string = manifest_file.read()
        manifest = json.loads(json_string)
    return manifest["component"], DatabaseRevision(manifest["revision"])


class CommecDatabaseUpdater:
    """
    Utility class, performs an update for a particular database.
    Handles fetching, writing, and version management.
    """

    def __init__(self, existing_location: os.PathLike, base_url: str):
        self.name = None
        self.write_location = existing_location
        self.base_url = base_url
        # 0.0 reads as "not installed", which is exactly what an unreadable
        # manifest tells us. Leaving this None makes every later .invalid()
        # call an AttributeError, which aborts a screen mid-setup.
        self.existing_revision = DatabaseRevision()
        try:
            self.name, self.existing_revision = read_manifest(existing_location)
        except json.JSONDecodeError as e:
            print(
                f"{C_F_ORANGE}{C_BOLD}X{C_RESET} Issue reading {existing_location}/manifest.json : {e}"
            )
        self.fetch_location = None
        self.update_required = True
        self.__update_message = None

        self.requested_revision = DatabaseRevision()

    @property
    def update_message(self):
        return (
            self.__update_message
            or f'{C_F_ORANGE}Unidentified{C_RESET} database "{self.name}" has no update option.'
        )

    def check_for_update(self, requested_revision: str):
        self.requested_revision = DatabaseRevision(requested_revision)

        # Database not supported.
        min_revision, max_revision = fetch_supported_revisions().get(
            self.name, (None, None)
        )
        if min_revision:
            if self.requested_revision >= max_revision:
                print(
                    f"{WARNING_CHECK}Version {VERSION} of commec supports "
                    f"to {self.name} <{str(max_revision)} and does not support"
                    f" revision {str(self.requested_revision)}.  "
                    f"\n   We recommend {C_F_ORANGE}you update commec{C_RESET}, and rerun this setup."
                )
                self.update_required = False
                self.__update_message = f"{C_F_ORANGE}commec package out of date{C_RESET} for latest {self.name} (revision {requested_revision})"
                return

        # Database does not yet exist.
        if not self.existing_revision or self.existing_revision.invalid():
            self.update_required = True
            self.__update_message = f"{C_F_ORANGE}Download{C_RESET} {self.name} (revision {requested_revision})"
            return

        # Database requires updating to latest revision
        if self.existing_revision < self.requested_revision:
            self.__update_message = f"{self.name} {self.existing_revision} will be {C_F_ORANGE}updated{C_RESET} to revision {self.requested_revision}."
            self.update_required = True
            return

        # Database has requested a downgrade to a specific revision
        if self.existing_revision > self.requested_revision:
            self.__update_message = f"{self.name} {self.existing_revision} will be {C_F_ORANGE}reverted{C_RESET} to revision {self.requested_revision}."
            self.update_required = True
            return

        # Database is up to date.
        self.__update_message = f"{self.name}: revision {self.existing_revision}. {C_F_ORANGE}Up to date!{C_RESET}"
        self.update_required = False
        return

    def perform_update(self):
        if not self.update_required:
            return

        os.makedirs(self.write_location, exist_ok=True)

        # print(f"{C_F_GRAY} Removing existing {self.name} ...")
        # Get a list of existing files at the location

        print(
            f"{STEP}Downloading {self.name} revision {str(self.requested_revision)}... "
        )
        # Calculate fetch_location
        self.fetch_location = os.path.join(self.name, str(self.requested_revision), "")
        tar_fetch_location = os.path.join(self.fetch_location, f"{self.name}.tar.zst")

        manifest_write_location = os.path.join(self.write_location, "manifest.json")
        tar_write_location = os.path.join(self.write_location, f"{self.name}.tar.zst")

        # Download tar from R2 self.fetch_location.
        download_success = save_r2_object(
            tar_fetch_location, tar_write_location, base_url=self.base_url
        )

        if not download_success:
            print(f"{ERROR_CHECK}Failed to download {self.name} ... ")
            return

        # Extract the Database tar.
        # Stdlib tarfile only gained native zstd support ('r:zst') in
        # Python 3.14 (PEP 784); this project caps at 3.13 due snakemake
        # benchmarking compatibility so we shell out to a locally installed tar
        # binary (with zstd support) rather than pull in a Python version dependency.
        print(f"{STEP}Extracting {self.name} database ... ")
        output = subprocess.run(
            ["tar", "--zstd", "-xf", tar_write_location, "-C", self.write_location]
        )
        if output.returncode > 0:
            print(f"{ERROR_CHECK}Error during Extraction: {output.stdout}")
            return

        # Remove the tar
        # os.remove(tar_write_location)

        # Write updated local manifest.
        manifest_data = {
            "component": str(self.name),
            "revision": str(self.requested_revision),
        }
        with open(manifest_write_location, mode="w", encoding="utf-8") as manifest_json:
            json.dump(manifest_data, manifest_json, indent=2)

        return


def create_updaters_from_config(
    config: dict, base_url: str
) -> dict[str, CommecDatabaseUpdater]:
    """
    Given a yaml config dict, create a named dict of updaters.
    """
    updaters = {}
    for database_name, database_info in config.get("databases", "").items():
        # The biorisk database points directly to the .hmm, whereas other
        # databases don't. It is best to just detect suffix, and remove, as
        # checking on name would affect other databases too.
        fileless_path = remove_filename_from_path(database_info.get("path"))
        updaters[database_name] = CommecDatabaseUpdater(fileless_path, base_url)
        updaters[database_name].name = database_name
    return updaters


@dataclass
class DatabaseRevision:
    """
    We don't have full semantic verisioning for databases, instead we use a
    revision system, i.e. 1.2, where the first number increments on breaking
    changes that require a commec update, and minor revisions increment the
    second number.
    """

    major: int = 0
    minor: int = 0

    def __init__(self, input_string: str = str("0.0")):
        extraction = str(input_string).split(".")
        if len(extraction) == 2:
            self.major = int(extraction[0])
            self.minor = int(extraction[1])
            return
        raise ValueError('Expected "X.X" string as input for Database Revision.')

    def invalid(self):
        return self.major == 0 and self.minor == 0

    def __eq__(self, other):
        return self.major == other.major and self.minor == other.minor

    def __lt__(self, other):
        if not isinstance(other, DatabaseRevision):
            return NotImplemented
        if self.major == other.major:
            return self.minor < other.minor
        return self.major < other.major

    def __gt__(self, other):
        if not isinstance(other, DatabaseRevision):
            return NotImplemented
        if self.major == other.major:
            return self.minor > other.minor
        return self.major > other.major

    def __le__(self, other):
        if not isinstance(other, DatabaseRevision):
            return NotImplemented
        if self.major == other.major:
            return self.minor <= other.minor
        return self.major <= other.major

    def __ge__(self, other):
        if not isinstance(other, DatabaseRevision):
            return NotImplemented
        if self.major == other.major:
            return self.minor >= other.minor
        return self.major >= other.major

    def __str__(self):
        return f"{self.major}.{self.minor}"


def resolve_base_url(
    cli_base_url: str | None,
    config: dict,
    config_source: str = "base_url in the config",
) -> str:
    """
    Work out which host to fetch databases from: a --base-url given on the
    command line wins, otherwise the config's base_url. `config_source` names
    the config in the message if its value is unusable, since the config may be
    one the user did not name explicitly.

    Any config merged from the packaged defaults carries a base_url already; the
    fallback is for one assembled by hand, which has no defaults behind it.
    """
    if cli_base_url is not None:
        return YamlIO.validate_base_url(cli_base_url, "--base-url")
    if config.get("base_url") is None:
        return YamlIO.default_base_url()
    return YamlIO.validate_base_url(config["base_url"], config_source)


def normalized_base_url(value, source: str) -> str | None:
    """
    Normalize a base URL for comparison, or None if it isn't one we could use.

    For a value that is only being consulted, not obeyed -- whether a config
    already records the URL in play. Anything that must reject a bad value calls
    `validate_base_url` directly, so the user hears about it.
    """
    if value is None:
        return None
    try:
        return YamlIO.validate_base_url(value, source)
    except ValueError:
        return None


def unwritable_download_directories(
    updaters: dict[str, "CommecDatabaseUpdater"],
) -> list[str]:
    """
    Which of the directories about to be downloaded into cannot be written to.

    Checked per database rather than only for the parent directory, since a config
    may place one database on a separate volume. A download that only finds out at
    extraction time has already spent the transfer.

    Only databases actually being updated are checked: an install that is already
    up to date writes nothing, so a directory it never touches being read-only is
    not this run's problem.
    """
    unwritable = []
    checked = set()
    for updater in updaters.values():
        if not updater.update_required:
            continue
        directory = str(updater.write_location)
        if directory in checked:
            continue
        checked.add(directory)
        if not check_directory_is_writable(directory):
            unwritable.append(directory)
    return unwritable


def warn_base_url_not_recorded(
    config_filepath: os.PathLike | str, base_url: str, recorded_base_url: str | None
) -> None:
    """
    Say how to make a base URL stick, when the config recording this install does
    not record it and this run isn't the one writing that config.

    A later run given that config with -y takes the URL from it, so leaving the two
    disagreeing means updates fetch from a host the databases did not come from.

    `recorded_base_url` is what that config currently names, already normalized,
    so it compares like for like against the URL actually used. Re-reading the
    file here instead would compare a raw value against a normalized one and
    report a config that differs only by a trailing slash as not recording it.
    """
    if recorded_base_url == base_url:
        return

    consequence = (
        f"will fetch from {recorded_base_url} instead"
        if recorded_base_url
        else "will not be able to read a host from it"
    )
    # Neutral about whether the file already has a base_url line: once a config
    # is merged with the packaged defaults it always names a host, so there is
    # no longer any way to tell an explicit default from an inherited one.
    print(
        f"{WARNING_CHECK}{config_filepath} does not record the base URL used here."
        f"\n   Set base_url in it to the following, or a later run given it with"
        f"\n   -y, including automatic updates during 'commec screen', {consequence}:"
        f"\n     base_url: {base_url}"
    )


def connection_failure_notes(err: BaseException, url: str) -> str | None:
    """
    Return some notes to print alongside a failed download, or None if we have
    nothing useful to add.
    """
    # urllib wraps the underlying error as URLError.reason. An ssl.SSLError also
    # has a .reason, but it's a string ("CERTIFICATE_VERIFY_FAILED"), so the
    # exception itself has to be checked as well as what it wraps.
    if not isinstance(err, ssl.SSLError) and not isinstance(
        getattr(err, "reason", None), ssl.SSLError
    ):
        return None

    host = urlparse(url).hostname or ""
    return (
        f"\n{C_F_GRAY}   Notes: this is a TLS/certificate failure rather than an HTTP\n"
        "   error, so the connection did not get far enough to download anything.\n"
        "   Possible causes, roughly in order of how often they turn up:\n"
        "     - A proxy, firewall, or antivirus product on the network is inspecting\n"
        "       HTTPS traffic and presenting its own certificate. If so, its root\n"
        "       certificate needs to be trusted (see SSL_CERT_FILE below), or the\n"
        "       host needs to be exempted from inspection.\n"
        "     - The CA certificates available to Python are out of date or\n"
        "       incomplete, which conda environments are especially prone to.\n"
        "     - A genuine problem with the certificate being served.\n"
        "   To see which certificate is actually arriving:\n"
        f"     openssl s_client -connect {host}:443 -servername {host} </dev/null \\\n"
        "       | openssl x509 -noout -issuer -subject -dates\n"
        "   To point commec at a specific CA bundle:\n"
        "     SSL_CERT_FILE=/path/to/ca-bundle.pem commec setup ...\n"
        f"{C_RESET}"
    )


def fetch_r2_object(
    object_path: str,
    base_url: str,
    timeout: int = 10,
    max_retries: int = 3,
    retry_delay: float = 2.0,
) -> bytes | None:
    """
    CURRENTLY only used for the latest.json, we need to use something more robust for
    the larger files.

    Download a single object from the R2 bucket, retrying transient
    network failures up to `max_retries` times with a short backoff.

    `object_path` is the object's key/path within the bucket and is joined onto
    `base_url`, which callers must supply -- it comes from `resolve_base_url`,
    so it is either the public default or whatever host the user pointed
    `commec setup --base-url` at.

    Returns the raw bytes of the object, or None if every attempt failed.
    """
    url = f"{base_url.rstrip('/')}/{object_path.lstrip('/')}"
    req = request.Request(url, headers={"User-Agent": "commec-setup"})

    for attempt in range(1, max_retries + 1):
        try:
            with request.urlopen(req, timeout=timeout) as response:
                if response.status == 200:
                    return response.read()
                print(
                    f"{ERROR_CHECK}HTTP error fetching {url}: status {response.status}"
                )
        except error.HTTPError as e:
            print(f"{ERROR_CHECK}HTTP Error fetching {url}: {e.code} - {e.reason}")
            if 400 <= e.code < 500:
                break  # Client errors (404, 403, ...) won't succeed on retry.
        except error.URLError as e:
            print(f"{ERROR_CHECK}URL Error fetching {url}: {e.reason}")
            # Only once the retries are spent, so it isn't repeated per attempt.
            if attempt == max_retries:
                notes = connection_failure_notes(e, url)
                if notes:
                    print(notes)

        if attempt < max_retries:
            time.sleep(retry_delay * attempt)

    return None


def save_r2_object(
    object_path: str,
    destination_path: os.PathLike | str,
    base_url: str,
    max_retries: int = 3,
    retry_delay: float = 2.0,
    chunk_size: int = 8 * 1024 * 1024,
) -> bool:
    """
    Download a (potentially very large, >7 GB) object from the R2 bucket and
    stream it to disk at `destination_path`, rather than holding it in
    memory like `fetch_r2_object`.

    A SHA-256 checksum is expected to exist alongside the object, at
    `object_path` with ".sha256" appended, and is fetched via
    `fetch_r2_object`. The object is hashed as it's streamed to disk, and
    the result is compared against this checksum, with a mismatch
    triggering a retry. A missing checksum fails the download, rather than
    leaving an unverified file in place.

    A failed download, or a checksum mismatch, is retried up to
    `max_retries` times.

    Returns True if the file was downloaded and verified, False otherwise.
    """
    expected_sha256_raw = fetch_r2_object(f"{object_path}.sha256", base_url=base_url)
    expected_sha256 = None
    if expected_sha256_raw is None:
        print(f"{ERROR_CHECK}Could not fetch checksum for {object_path}")
        return False
    else:
        expected_sha256 = expected_sha256_raw.decode("utf-8").strip().split()[0]

    url = f"{base_url.rstrip('/')}/{object_path.lstrip('/')}"
    req = request.Request(url, headers={"User-Agent": "commec-setup"})
    destination_path = Path(destination_path)

    for attempt in range(1, max_retries + 1):
        try:
            with request.urlopen(req, timeout=None) as response:
                if response.status != 200:
                    print(
                        f"{ERROR_CHECK}HTTP error fetching {url}: status {response.status}"
                    )
                else:
                    total_size = int(response.getheader("Content-Length") or 0)
                    downloaded = 0
                    last_progress_time = 0.0
                    progress_interval = 1 / 20  # throttle to at most ~20 updates/sec

                    sha256_hash = hashlib.sha256()
                    with open(destination_path, "wb") as out_file:
                        while chunk := response.read(chunk_size):
                            out_file.write(chunk)
                            sha256_hash.update(chunk)
                            downloaded += len(chunk)

                            now = time.monotonic()
                            if now - last_progress_time >= progress_interval:
                                print_progress(downloaded, total_size)
                                last_progress_time = now

                    # Always print a finished bar!
                    print_progress(downloaded, total_size)
                    print()

                    if sha256_hash.hexdigest() == expected_sha256:
                        return True
                    print(
                        f"{ERROR_CHECK}Checksum mismatch for {url}: "
                        f"expected {expected_sha256}, got {sha256_hash.hexdigest()}"
                    )
        except error.HTTPError as e:
            print(f"{ERROR_CHECK}HTTP Error fetching {url}: {e.code} - {e.reason}")
            if 400 <= e.code < 500:
                break  # Client errors (404, 403, ...) won't succeed on retry.
        except (error.URLError, OSError) as e:
            print(f"{ERROR_CHECK}Error fetching {url}: {e}")
            if attempt == max_retries:
                notes = connection_failure_notes(e, url)
                if notes:
                    print(notes)

        if attempt < max_retries:
            time.sleep(retry_delay * attempt)

    destination_path.unlink(missing_ok=True)
    return False


def print_progress(downloaded: int, total_size: int) -> None:
    """
    Render an in-place, single-line progress indicator to stdout.

    Clears only the current terminal line (rather than scrolling) before
    redrawing, so repeated calls animate a single line in place.

    `total_size` of 0 means the size couldn't be determined (e.g. no
    Content-Length header on the response) -- a percentage/bar can't be
    computed without a denominator, so the raw amount downloaded is
    reported instead, to keep terminal output compact.
    """
    downloaded_mb = downloaded / (1024 * 1024)

    if total_size <= 0:
        line = f"{BULLET}Downloading: {downloaded_mb:.1f} MB"
    else:
        total_mb = total_size / (1024 * 1024)
        fraction = min(downloaded / total_size, 1.0)
        suffix = f" {fraction * 100:5.1f}% ({downloaded_mb:7.1f}/{total_mb:7.1f} MB)"
        bar_width = max(80 - len(suffix) - 2, 10)
        filled = int(bar_width * fraction)
        bar = (
            C_F_ORANGE + "●" * filled + C_F_BLUE + "•" * (bar_width - filled) + C_RESET
        )
        line = f"[{bar}]{C_F_GRAY}{suffix}{C_RESET}"

    sys.stdout.write(f"\033[2K\r{line}")
    sys.stdout.flush()


def read_config_yaml_for_database_setup(
    config_yaml_filepath: os.PathLike | str, database_dir: os.PathLike | str = None
):
    """
    Reads a config yaml, updated from the defaults, and parses the output for
    the control_list directory as per a commec screen run. Used instead of passing
    the directory of the control list

    :param config_yaml_filepath: Description
    :type config_yaml_filepath: os.PathLike | str
    :param database_dir: The parent directory named by `-d`, when one was given.
        Wins over the config's own base_paths.default, so a config may leave that
        unset and be shared across machines that keep databases in different
        places. Paths the config states outright are left alone, which is how a
        database is kept on a separate volume.
    """

    output_config = None

    # Read package-level configuration defaults
    default_yaml = importlib.resources.files("commec").joinpath(
        DEFAULT_CONFIG_YAML_PATH
    )
    if default_yaml.exists():
        output_config = YamlIO.load_config_from_yaml(str(default_yaml))
    else:
        raise FileNotFoundError(
            f"No default yaml found. Expected at {DEFAULT_CONFIG_YAML_PATH}"
        )

    # Override configuration with any in user-provided YAML file. A path that
    # isn't there is named as such, rather than left to surface later as a
    # puzzling complaint about a key the config never had a chance to set.
    # `commec screen` reports the same mistake the same way, and as the same kind
    # of error, so the CLI turns both into one line instead of a traceback.
    if not os.path.exists(config_yaml_filepath):
        raise YamlIOValidationError(f"--config YAML not found: {config_yaml_filepath}")

    output_config = YamlIO.update_config_from_yaml(output_config, config_yaml_filepath)

    if database_dir is not None:
        output_config["base_paths"]["default"] = str(Path(database_dir).resolve())

    output_config = YamlIO.format_config_paths(output_config)

    return output_config


def check_directory_is_writable(input_directory: str) -> str:
    """
    Checks a directory is viable by
    * Expanding terminal variables, user, and resolving the full path and
    * Checking if it exists or
    * Creating it and destroying it.

    It returns a str representing the both the truthiness of the outcome,
    as well as the valid path.
    """
    path = Path(os.path.expandvars(input_directory))

    # Catches accidental ~/commec-dbs/ vs ~commec-dbs/ when theres no commec-dbs username.
    try:
        path = path.expanduser()
    except RuntimeError:
        print(
            "User expansion for path failed, ensure you are using"
            ' "~/" for self, or a valid user with "~username/".'
        )
        return ""

    try:
        path = path.resolve()
    except RuntimeError:
        return ""

    if path.exists():
        return path

    if path.is_reserved():
        print("This path contains reserved characters for this Operating System.")
        return ""

    # Handily, all sorts of special characters are identified with a %XX, within posix, and are replaced
    # by similar characters during mkdir, whilst technically legal, lets recommend against cursed dir names.
    if "%" in path.as_posix():
        print(
            'Please avoid using special characters ("|}{":?><*&" etc) in filepath names.'
        )
        return ""

    # If the path doesn't exist, the best way to know if user input is valid, is to try make it.
    # Find the part of the directory which is new, so we can delete only it after.
    path_to_remove_dirs = Path(path.parts[0])
    for part in path.parts:
        if path_to_remove_dirs.exists():
            path_to_remove_dirs = path_to_remove_dirs / part
            continue
        break

    # Create the directory, and delete anything created. Failing to create it is
    # the answer this function exists to give, so it is reported by returning
    # falsy like every other rejection here, rather than as a raised OSError that
    # callers checking the return value would never think to catch.
    try:
        os.makedirs(path, exist_ok=True)
    except OSError:
        return ""
    if path.exists():
        try:
            shutil.rmtree(path_to_remove_dirs)
        except OSError:
            pass
        return path
    return ""


def check_for_updates(config: dict) -> tuple[bool, dict[str, CommecDatabaseUpdater]]:
    """
    Helper function for calling from commec, quickly assess if an update is
    required for yaml config dict, and output the updaters of those databases
    which require updates. Call perform_update() on updates to complete updates.
    Requests latest updates only.

    Fetches from the config's base_url if it records one (as written by
    `commec setup --base-url`), so an update triggered from a screen given that
    config with -y uses the same host the databases were installed from. A screen
    given only -d has no config to read a host from, and falls back to the default;
    automatic updates are off unless a config turns them on, so the two go together.
    """
    base_url = resolve_base_url(None, config)
    updaters = create_updaters_from_config(config, base_url)
    revisions = fetch_latest_revisions(base_url)
    if not revisions:
        # Unreachable host, or a latest.json that wouldn't parse. Nothing to
        # compare against, so leave the installed databases as they are.
        return False, {}
    latest_revisions = revisions["latest"]
    for dbname, latest in latest_revisions.items():
        db_to_update = updaters.get(dbname)
        if not db_to_update:
            continue  # Ignore yaml databases that aren't in latest revisions.
        db_to_update.name = dbname
        db_to_update.check_for_update(latest)

    update_flag: bool = any(
        [
            (db.update_required and not db.existing_revision.invalid())
            for db in updaters.values()
        ]
    )
    return update_flag, updaters


def run(arguments):
    """Entry point for CommecSetup"""
    print(SPLASH_IMAGE)
    print(SPLASH_TEXT)
    CommecSetup(arguments)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    args = parser.parse_args()
    run(args)
