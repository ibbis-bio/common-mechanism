#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Abstract base class defining a shared interface for search tools.
"""
from abc import ABC, abstractmethod
import os
from dataclasses import dataclass
import subprocess
import logging

logger = logging.getLogger(__name__)

@dataclass
class SearchToolVersion:
    """Container class for outputting version related information from a database."""

    tool_info: str = "x.x.x"
    database_info: str = "x.x.x"


class DatabaseValidationError(Exception):
    """Custom exception for database validation errors."""

class SearchHandler(ABC):
    """
    Abstract class defining tool interface including a database directory / file to search, an input
    query, and an output file to be used for screening.
    """

    def __init__(
        self,
        database_file: str | os.PathLike,
        input_file: str | os.PathLike,
        out_file: str | os.PathLike,
        **kwargs,
    ):
        """
        Initialise a Search Handler.

        Parameters
        ----------
        database_file : str | os.PathLike
            Path to the database file.
        input_file : str | os.PathLike
            Path to the input file to be processed.
        out_file : str | os.PathLike
            Path where the output will be saved.

        Keyword Arguments
        -----------------
        threads : int, optional
            Number of threads to use for processing. Default is 1.
        force : bool, optional
            Whether to force overwrite existing files. Default is False.

        Notes
        -----
        - `database_file`, `input_file`, and `out_file` are validated on instantiation.
        """

        self.db_file = os.path.abspath(os.path.expanduser(database_file))
        self.input_file = os.path.abspath(os.path.expanduser(input_file))
        self.out_file = os.path.abspath(os.path.expanduser(out_file))
        self.threads = kwargs.get('threads', 1)
        self.force = kwargs.get('force', False)
        self.arguments_dictionary = {}
        self.successful = True

        # Only validate database files if we actually intend on using them
        if not self.should_use_existing_output:
            self._validate_db()

        self.version_info = self.get_version_information()

    @property
    def db_directory(self):
        """Directory where databases to be searched are located."""
        return os.path.dirname(self.db_file)

    @property
    def temp_log_file(self):
        """Temporary log file used for this search. Based on outfile name."""
        return f"{self.out_file}.log.tmp"

    @property
    def should_use_existing_output(self) -> bool:
        """
        True if (1) search is not forced and (2) output exists and is valid.
        """
        return not self.force and self.validate_output()

    def search(self):
        """
        Wrapper for _search, skipping if existing output should not be overwritten.
        """
        if self.should_use_existing_output:
            logger.warning("%s expected output data already exists, "
                         "will use existing data found in:",
                         self.__class__.__name__)
            logger.warning(self.out_file, extra = {"no_prefix" : True, "cap":True})
        else:
            self._search()

    @abstractmethod
    def _search(self):
        """
        Use a tool to search the input query against a database.
        Should be implemented by all subclasses to perform the actual search against the database.
        """

    @abstractmethod
    def read_output(self):
        """
        Returns the output of the handler in the form of a pandas dataframe.
        """

    @abstractmethod
    def get_version_information(self) -> SearchToolVersion:
        """
        Provide version for the search tool used, to allow reproducibility.
        This method should be implemented by all subclasses to return tool-specific version info.
        """

    def validate_output(self):
        """
        Check the output file contains something, indicating that the search ran.
        Can be overridden if more complex checks for a particular tool are desired.
        Is overridden for Diamond outputs, which have no header information, and simply only
        checks for file-existance, rather than lack of content, for example.
        """
        return not self.has_empty_output()

    def _validate_db(self):
        """
        Validates that the database directory and file exists. Called on init.
        """
        if not os.path.isdir(self.db_directory):
            raise DatabaseValidationError(
                f"Screening database directory not found at: {self.db_directory}."
                " Screening directory path can be set via --databases option or --config yaml."
            )

        if not os.path.isfile(self.db_file):
            raise DatabaseValidationError(
                f"Provided database file not found: {self.db_file}."
                " File location can be set via --databases option or --config yaml."
            )

    def has_empty_output(self) -> bool:
        """Check if the output file is empty or non-existent."""
        try:
            return os.path.getsize(self.out_file) == 0
        except OSError:
            # Errors such as FileNotFoundError considered empty
            return True

    def has_hits(self) -> bool:
        """Check if the output file has any hits (lines that do not start with '#')."""
        try:
            with open(self.out_file, "r", encoding="utf-8") as file:
                return any(not line.strip().startswith("#") for line in file)
        except FileNotFoundError:
            return False

    def format_args_for_cli(self) -> list:
        """
        Format `self.arguments_dictionary` into a list of strings for use in the command line.
        """
        formatted_args = []
        for key, value in self.arguments_dictionary.items():
            formatted_args.append(str(key))
            if isinstance(value, list):
                formatted_args.append(" ".join(map(str, value)))
            elif value is not None:
                formatted_args.append(str(value))
        return formatted_args

    def run_as_subprocess(self, command, out_file, raise_errors=False):
        """
        Run a command using subprocess.run, piping stdout and stderr to `out_file`.
        """
        self.successful = False

        logger.debug("SUBPROCESS: %s", " ".join(command))
        logger.debug(" ".join(command), extra = {"no_prefix":True,"cap":True})

        with open(out_file, "a", encoding="utf-8") as f:
            result = subprocess.run(
                command, stdout=f, stderr=subprocess.STDOUT, check=raise_errors
            )

            if result.returncode != 0:
                command_str = " ".join(command)
                logger.error(
                    "\t command '%s' failed with error '%s'",
                    command_str,
                    result.stderr,
                )
                raise RuntimeError(
                    f"subprocess.run of command '{command_str}' encountered error."
                    f" Check {out_file} for logs."
                )
            
        self.successful = True

    def __del__(self):
        if os.path.exists(self.temp_log_file) and self.successful:
            os.remove(self.temp_log_file)