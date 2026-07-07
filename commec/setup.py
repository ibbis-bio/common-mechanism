#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Module for CLI setup of Commec, such that required 
databases are downloaded in a desired database directory.
"""
import sys
import os
import shutil
import argparse
import subprocess
import time
import hashlib
from pathlib import Path
import ftplib
import importlib.resources
from urllib import request, error, parse
import zipfile
import yaml
import json
from yaml.parser import ParserError
from dataclasses import dataclass

from commec.config.constants import DEFAULT_CONFIG_YAML_PATH
from commec.config import yaml_io as YamlIO
from commec.utils.file_utils import expand_and_normalize

from  commec import __version__ as VERSION
DESCRIPTION = """Helper script for downloading or updating the databases
 required for running the Common Mechanism Screen"""

C_F_ORANGE = "\033[38;5;202m" # Colour Foreground Orange.
C_F_GRAY = "\033[38;5;242m" # Colour Foreground Gray
C_F_BLUE = "\033[38;5;17m" # Colour Foreground Blue
C_B_BLUE = "\033[48;5;17m" # Colour Background Blue
C_RESET = "\033[0m" # Reset Console Formatting.
C_BOLD = "\033[1m" # Reset Console Formatting.

ERROR_CHECK = C_F_ORANGE + C_BOLD + " X " + C_RESET
STEP = " ➔ "
BULLET = "  ● "

# Base URL of the R2 bucket that hosts the commec databases and their
# version manifest. Currently the public dev endpoint; swap for a stable
# custom domain once one is provisioned (r2.dev is not meant for
# production traffic -- no SLA, unpublished rate limits).
R2_PUBLIC_BASE_URL = "https://pub-411843ff50b843c6b1e8d9600eed9093.r2.dev"

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

    # These two should be a mutually exclusive group.
    exclusive_group = input_options.add_mutually_exclusive_group()
    exclusive_group.add_argument(
        "-d",
        "--databases",
        dest="database_dir",
        default="",
        help="Path to a parent directory to download, or update, the Commec"
        " databases.",
    )
    exclusive_group.add_argument(
        "-y",
        "--config",
        dest="yaml_file",
        default="",
        help="Alternative to -d, --databases. Specify databases location with a commec yaml configuration file",
    )

    # Optional
    input_options.add_argument(
        "-r",
        "--dryrun",
        dest="dry_run",
        default=False,
        action="store_true",
        help="Output potential update information, without performing any download or update.",
    )

    return input_options

class CommecSetup:
    def __init__(self, args):
        working_directory = Path()
        # Updaters dictionary, of {[name : UpdaterObject]}
        updaters = {}
        self.yaml_mode = False
        self.config = YamlIO.get_defaults()

        # Derive databases directory / yaml outputs.abs
        if args.database_dir:
            working_directory = expand_and_normalize(args.database_dir)
        if args.yaml_file:
            self.config = read_config_yaml_for_database_setup(args.yaml_file)
            working_directory = expand_and_normalize(self.config.get("base_paths").get("default"))
            self.yaml_mode = True
            for database_name, database_info in self.config.get("databases","").items():
                # The biorisk database points directly to the .hmm, whereas other
                # databases don't. It is best to just detect suffix, and remove, as 
                # checking on name would affect other databases too.
                db_path = Path(database_info.get("path"))
                fileless_path = db_path.parent if db_path.suffix else db_path
                updaters[database_name] = CommecDatabaseUpdater(fileless_path)
                updaters[database_name].name = database_name

        # Ensure working_directory is a valid path (even if it doesn't exist yet)
        if not check_directory_is_writable(working_directory):
            print(f"{ERROR_CHECK}Failed to parse working directory: {working_directory}")

        # Fetch latest.json from R2.abs
        print(f"{STEP}Fetching latest commec database revision information ...")
        latest = fetch_latest_revisions()

        if not latest:
            print(f"{ERROR_CHECK} Could not fetch latest revisions. Check internet connection and try again.")
            return

        # For each database, create an updater.
        databases = latest["databases"]
        for database_name, revision in databases.items():
            updater = updaters.get(database_name)
            if not updater:
                download_location = os.path.join(working_directory,database_name)
                updater = CommecDatabaseUpdater(download_location)
                updaters[database_name] = updater
                if self.yaml_mode:
                    print(f"{ERROR_CHECK}{database_name} not present in provided configuration yaml, it may be out of date.")
            updater.name = database_name
            updater.check_for_update(revision)

        # Check or ophaned databases:
        for db_name, updater in updaters.items():
            if not db_name in databases.keys():
                print(f"{ERROR_CHECK}{db_name} has no formal update route. Input yaml has deprecated database entries.")
                updater.update_required = False

        # Each updater, reports on its status, does it need to update?
        print(f"{STEP}The following actions have been identified ...")
        for _db_name, updater in updaters.items():
            print(BULLET,updater.update_message)

        # Log output of intended changes.
        if args.dry_run:
            print(
                f"\n{C_F_GRAY}This was a dry run. Run {C_F_ORANGE}\'commec"
                f" setup\'{C_F_GRAY} again without {C_F_ORANGE}\'-r\'{C_F_GRAY}"
                f" to continue with the update.{C_RESET}")
            return

        # Perform the database updates, consider making this async.
        os.makedirs(working_directory, exist_ok=True)
        for database_name, updater in updaters.items():
            updater.perform_update()

        # Create example config.yaml if we only passed a directory.
        if not self.yaml_mode:
            self.config["base_paths"]["default"] = str(Path(working_directory).resolve())
            output_yaml_filepath = os.path.join(working_directory, "config.yaml")
            with open(output_yaml_filepath, 'w', encoding="utf-8") as output_config_file:
                yaml.safe_dump(self.config, output_config_file)
            print(f"{STEP}An example config.yaml for commec was created in the provided directory.")

        print(f"{STEP}Update check complete! Have a Biosafe-and-secure day!")

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
    data_filename = importlib.resources.files("commec").joinpath("database_compatibility.yaml")
    supported_revision_data = YamlIO.load_config_from_yaml(data_filename)
    for db_name, revision_string in supported_revision_data["supported_database_revisions"].items():
        # Expects the following str format: ">=1.0,<2.0"
        major_min = revision_string.split(",")[0].split(".")[0][2:]+".0"
        major_max = revision_string.split(",")[1].split(".")[0][1:]+".0"
        supported_revision_data[db_name] = (DatabaseRevision(major_min), DatabaseRevision(major_max))
    return supported_revision_data

def fetch_latest_revisions() -> dict | None:
    """
    A latest.json manifest exists at the root of the R2 bucket, listing
    the latest revision of each database. Downloads and parses it into a
    dict. Returns None if it could not be fetched or parsed.
    """
    raw = fetch_r2_object("latest.json")
    if raw is None:
        return None

    try:
        return json.loads(raw.decode("utf-8"))
    except (json.JSONDecodeError, UnicodeDecodeError) as e:
        print(f"Invalid JSON in latest.json manifest: {e}")
        return None

class CommecDatabaseUpdater:
    """
    Utility class, performs an update for a particular database.
    Handles fetching, writing, and version management.
    """
    def __init__(self, existing_location : os.PathLike):
        self.name = None
        self.write_location = existing_location
        self.existing_revision = None

        self.read_manifest()

        self.fetch_location = None
        self.update_required = True
        self.__update_message = None

        self.latest_revision = DatabaseRevision()
        self.md5 = ""
        self.sha256 = ""
        self.size_bytes = ""

    @property
    def update_message(self):
        return (self.__update_message or
                f"{C_F_ORANGE}Unidentified{C_RESET} database \"{self.name}\" has no update option.")

    def read_manifest(self):
        """
        All databases should have a manifest.json file.
        """
        manifest_filename = os.path.join(self.write_location, "manifest.json")
        if os.path.exists(self.write_location) and os.path.isfile(manifest_filename):
            try:
                with open(manifest_filename, "r", encoding="utf-8") as manifest_file:
                    json_string = manifest_file.read()
                manifest = json.loads(json_string)
                self.name = manifest["component"]
                self.existing_revision = DatabaseRevision(manifest["revision"])
            except json.JSONDecodeError as e:
                print(f"{C_F_ORANGE}{C_BOLD}X{C_RESET} Issue reading {manifest_filename} : {e}")
        return

    def check_for_update(self, latest_revision : str):
        self.latest_revision = DatabaseRevision(latest_revision)

        # Database not supported.
        min_revision, max_revision = fetch_supported_revisions().get(self.name, (None,None))
        if min_revision:
            if self.latest_revision >= max_revision:
                print(
                    f"{ERROR_CHECK}Version {VERSION} of commec supports "
                    f"to {self.name} <{str(max_revision)} and does not support"
                    f" revision {str(max_revision)}.  "
                    f"\n   We recommend {C_F_ORANGE}you update commec{C_RESET}, and rerun this setup.")
                self.update_required = False
                self.__update_message = f"{C_F_ORANGE}commec package out of date{C_RESET} for latest {self.name} (revision {latest_revision})"
                return

        # Database does not yet exist.
        if not self.existing_revision:
            self.update_required = True
            self.__update_message = f"{C_F_ORANGE}Download{C_RESET} {self.name} (revision {latest_revision})"
            return

        # Database requires updating to latest revision
        if self.existing_revision < self.latest_revision:
            self.__update_message = f"{self.name} ({self.existing_revision}) will be {C_F_ORANGE}updated{C_RESET} to revision {latest_revision}."
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

        #print(f"{C_F_GRAY} Removing existing {self.name} ...")
        # Get a list of existing files at the location

        print(f"{STEP}Downloading latest {self.name} ... ")
        # Calculate fetch_location
        self.fetch_location = os.path.join(self.name, str(self.latest_revision), "")
        manifest_fetch_location = os.path.join(self.fetch_location, "manifest.json")
        tar_fetch_location = os.path.join(self.fetch_location, f"{self.name}.tar.zst")

        manifest_write_location = os.path.join(self.write_location, "manifest.json")
        tar_write_location = os.path.join(self.write_location, f"{self.name}.tar.zst")

        # Download tar from R2 self.fetch_location.
        download_success = save_r2_object(tar_fetch_location, tar_write_location)
        
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
        #os.remove(tar_write_location)

        # Write updated local manifest.
        manifest_data = {
            "component" : str(self.name),
            "revision" : str(self.latest_revision),
        }
        with open(manifest_write_location, mode='w', encoding = "utf-8") as manifest_json:
            json.dump(manifest_data, manifest_json, indent = 2)

        return

@dataclass
class DatabaseRevision:
    """
    We don't have full semantic verisioning for databases, instead we use a 
    revision system, i.e. 1.2, where the first number increments on breaking
    changes that require a commec update, and minor revisions increment the
    second number.
    """
    major : int = 0
    minor : int = 0
    def __init__(self, input_string : str = str("0.0")):
        extraction = input_string.split(".")
        if len(extraction) == 2:
            self.major, self.minor = extraction
            return
        raise ValueError("Expected \"X.X\" string as input for Database Revision.")
    def __eq__(self, other):
        return (self.major == other.major 
                and self.minor == other.minor)
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

def fetch_r2_object(
    object_path: str,
    base_url: str = R2_PUBLIC_BASE_URL,
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
    `base_url`. `base_url` defaults to the public r2.dev endpoint but can
    be swapped for a stable custom domain fronting the same bucket without
    changing any calling code.

    Returns the raw bytes of the object, or None if every attempt failed.
    """
    url = f"{base_url.rstrip('/')}/{object_path.lstrip('/')}"
    req = request.Request(url, headers={"User-Agent": "commec-setup"})

    for attempt in range(1, max_retries + 1):
        try:
            with request.urlopen(req, timeout=timeout) as response:
                if response.status == 200:
                    return response.read()
                print(f"{ERROR_CHECK}HTTP error fetching {url}: status {response.status}")
        except error.HTTPError as e:
            print(f"{ERROR_CHECK}HTTP Error fetching {url}: {e.code} - {e.reason}")
            if 400 <= e.code < 500:
                break  # Client errors (404, 403, ...) won't succeed on retry.
        except error.URLError as e:
            print(f"{ERROR_CHECK}URL Error fetching {url}: {e.reason}")

        if attempt < max_retries:
            time.sleep(retry_delay * attempt)

    return None

def save_r2_object(
    object_path: str,
    destination_path: os.PathLike | str,
    base_url: str = R2_PUBLIC_BASE_URL,
    max_retries: int = 3,
    retry_delay: float = 2.0,
    chunk_size: int = 8 * 1024 * 1024,
) -> bool:
    """
    Download a (potentially very large, >7 GB) object from the R2 bucket and
    stream it to disk at `destination_path`, rather than holding it in
    memory like `fetch_r2_object`.

    An MD5 checksum is expected to exist alongside the object, at
    `object_path` with ".md5" appended, and is fetched via
    `fetch_r2_object`. If found, the object is hashed as it's streamed to
    disk, and the result is compared against this checksum, with a
    mismatch triggering a retry. If no checksum is found, the download is
    still attempted, just without any verification.

    A failed download, or a checksum mismatch, is retried up to
    `max_retries` times.

    Returns True if the file was downloaded successfully (and verified,
    if a checksum was available), False otherwise.
    """
    expected_sha256_raw = fetch_r2_object(f"{object_path}.sha256")
    expected_sha256 = None
    if expected_sha256_raw is None:
        print(f"{ERROR_CHECK}Could not fetch checksum for {object_path}, proceeding without verification.")
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
                    print(f"{ERROR_CHECK}HTTP error fetching {url}: status {response.status}")
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
        bar = C_F_ORANGE + "●" * filled + C_F_BLUE + "•" * (bar_width - filled) + C_RESET
        line = f"[{bar}]{C_F_GRAY}{suffix}{C_RESET}"

    sys.stdout.write(f"\033[2K\r{line}")
    sys.stdout.flush()

def read_config_yaml_for_database_setup(config_yaml_filepath : os.PathLike | str):
    """
    Reads a config yaml, updated from the defaults, and parses the output for
    the control_list directory as per a commec screen run. Used instead of passing
    the directory of the control list
    
    :param config_yaml_filepath: Description
    :type config_yaml_filepath: os.PathLike | str
    """

    output_config = None

    # Read package-level configuration defaults
    default_yaml = importlib.resources.files("commec").joinpath(DEFAULT_CONFIG_YAML_PATH)
    if default_yaml.exists():
        output_config = YamlIO.load_config_from_yaml(str(default_yaml))
    else:
        raise FileNotFoundError(
            f"No default yaml found. Expected at {DEFAULT_CONFIG_YAML_PATH}"
            )

    # Override configuration with any in user-provided YAML file
    if os.path.exists(config_yaml_filepath):
        output_config = YamlIO.update_config_from_yaml(output_config, config_yaml_filepath)

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
    except(RuntimeError):
        print("User expansion for path failed, ensure you are using"
                " \"~/\" for self, or a valid user with \"~username/\".")
        return ""

    try:
        path = path.resolve()
    except(RuntimeError):
        return ""

    if path.exists():
        return path
    
    if path.is_reserved():
        print("This path contains reserved characters for this Operating System.")
        return ""
    
    # Handily, all sorts of special characters are identified with a %XX, within posix, and are replaced
    # by similar characters during mkdir, whilst technically legal, lets recommend against cursed dir names.
    if '%' in path.as_posix():
        print("Please avoid using special characters (\"|}{\":?><*&\" etc) in filepath names.")
        return ""

    # If the path doesn't exist, the best way to know if user input is valid, is to try make it.
    # Find the part of the directory which is new, so we can delete only it after.
    path_to_remove_dirs = Path(path.parts[0])
    for part in path.parts:
        if path_to_remove_dirs.exists():
            path_to_remove_dirs = path_to_remove_dirs / part
            continue
        break

    # Create the directory, and delete anything created.
    os.makedirs(path, exist_ok=True)
    if path.exists():
        try:
            shutil.rmtree(path_to_remove_dirs)
        except OSError:
            pass
        return path
    return ""

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
