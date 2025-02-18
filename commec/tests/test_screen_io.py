import pytest
from unittest.mock import mock_open, patch
import os

from commec.config.io_parameters import ScreenIOParameters, IoValidationError
from commec.screen import add_args, ScreenArgumentParser


@pytest.fixture
def test_data_dir():
    return os.path.join(os.path.dirname(__file__), "test_data")


@pytest.fixture
def database_dir():
    return os.path.join(os.path.dirname(__file__), "test_dbs")


@pytest.mark.parametrize(
    "fasta_name",
    [
        "single_record.fasta",
        "multiple_records.fasta",
        "has_empty_record.fasta",
        "has_empty_description.fasta",
        "has_records_with_same_description.fasta",
    ],
)
def test_default_parameters(fasta_name, test_data_dir, database_dir, tmp_path):
    input_fasta = os.path.join(test_data_dir, fasta_name)
    with patch(
        "sys.argv",
        ["test.py", "-f", input_fasta, "-d", database_dir, "-o", str(tmp_path)],
    ):
        parser = ScreenArgumentParser()
        add_args(parser)
        screen_io = ScreenIOParameters(parser.parse_args())
        assert screen_io.setup()


@pytest.mark.parametrize(
    "fasta_name,expected_record_count",
    [
        pytest.param("single_record.fasta", 1),
        pytest.param("multiple_records.fasta", 2),
    ],
)
def test_parse_input_fasta(
    fasta_name, expected_record_count, test_data_dir, database_dir, tmp_path
):
    input_fasta = os.path.join(test_data_dir, fasta_name)
    with patch(
        "sys.argv",
        ["test.py", "-f", input_fasta, "-d", database_dir, "-o", str(tmp_path)],
    ):
        parser = ScreenArgumentParser()
        add_args(parser)
        screen_io = ScreenIOParameters(parser.parse_args())
        screen_io.setup()

    queries = screen_io.parse_input_fasta()
    assert len(queries) == expected_record_count


@pytest.mark.parametrize(
    "fasta_name",
    [
        "has_empty_record.fasta",
        "has_empty_description.fasta",
        "has_records_with_same_description.fasta",
    ],
)
def test_parse_invalid_input_fasta(fasta_name, test_data_dir, database_dir, tmp_path):
    input_fasta = os.path.join(test_data_dir, fasta_name)
    with patch(
        "sys.argv",
        ["test.py", "-f", input_fasta, "-d", database_dir, "-o", str(tmp_path)],
    ):
        parser = ScreenArgumentParser()
        add_args(parser)
        screen_io = ScreenIOParameters(parser.parse_args())
        screen_io.setup()

    with pytest.raises(IoValidationError):
        screen_io.parse_input_fasta()
