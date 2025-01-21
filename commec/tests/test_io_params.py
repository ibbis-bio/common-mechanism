import pytest
from unittest.mock import patch
import os
import argparse

from commec.config.io_parameters import ScreenIO
from commec.screen import add_args

INPUT_QUERY = os.path.join(os.path.dirname(__file__), "test_data/single_record.fasta")
DATABASE_DIRECTORY = os.path.join(os.path.dirname(__file__), "test_dbs/")

@patch("sys.argv", ["test_io_params.py", "-f", INPUT_QUERY, "-d", DATABASE_DIRECTORY])
def test_default_parameters():
    args = argparse.ArgumentParser()
    add_args(args)
    args = args.parse_args()
    new_params = ScreenIO(args)
    assert new_params.setup()
