#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Static functions useful for dealing with common file parsing tasks.
"""

import argparse
import os

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
