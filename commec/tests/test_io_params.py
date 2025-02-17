import pytest
from unittest.mock import patch
import os
import argparse

from commec.config.io_parameters import ScreenIOParameters
from commec.cli import ScreenArgumentParser
from commec.screen import add_args

INPUT_QUERY = os.path.join(os.path.dirname(__file__), "test_data/single_record.fasta")
DATABASE_DIRECTORY = os.path.join(os.path.dirname(__file__), "test_dbs/")

@patch("sys.argv", ["test_io_params.py", "-f", INPUT_QUERY, "-d", DATABASE_DIRECTORY, "-o", "output_test"])
def test_default_parameters(tmp_path):
    args = ScreenArgumentParser()
    add_args(args)
    args = args.parse_args()
    args.output_prefix = tmp_path
    new_params = ScreenIOParameters(args)
    assert new_params.setup()

