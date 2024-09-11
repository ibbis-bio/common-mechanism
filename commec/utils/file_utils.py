#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Module for FileTools, containing Static functions useful for dealing with common file parsing tasks.
"""

import argparse
import os

class FileTools():
    """
    Static function container, useful for dealing with common file parsing tasks.
    """

    # Can be moved somewhere else.
    @staticmethod
    def is_empty(filepath: str) -> bool:
        """
        is_empty

        usage: check that a file is empty
        input:
        - name of file
        """
        try:
            abspath = os.path.abspath(os.path.expanduser(filepath))
            return os.path.getsize(abspath) == 0
        except OSError:
            # If there is an error (including FileNotFoundError) consider it empty
            return True

    # Note this is typically used with a Blast outputfile, and thus may be better put in BlastTools.
    # Moved away to all databases.
    @staticmethod
    def has_hits(filepath: str) -> bool:
        """
        has_hits
        usage: check to see if the file contains any hits (lines that don't start with #)
        input:
        - path to file
        """
        try:
            with open(filepath, "r", encoding="utf-8") as file:
                for line in file:
                    # Strip leading and trailing whitespace and check the first character
                    if not line.strip().startswith("#"):
                        # Found a hit!
                        return True
            return False
        except FileNotFoundError:
            # The file does not exist
            return False
        


    # Below go to config parameters.
    @staticmethod
    def directory_arg(path):
        """Raise ArgumentTypeError if `path` is not a directory."""
        if not os.path.isdir(path):
            raise argparse.ArgumentTypeError(f"{path} is not a valid directory path")
        return path

    @staticmethod
    def file_arg(path):
        """Raise ArgumentTypeError if `path` is not a file."""
        if not os.path.isfile(path):
            raise argparse.ArgumentTypeError(f"{path} is not a valid file")
        if not os.path.getsize(path) > 0:
            raise argparse.ArgumentTypeError(f"{path} is an empty file")
        return path


    
    