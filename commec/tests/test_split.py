import os
from unittest.mock import mock_open, patch

import pytest
from Bio import SeqIO

from commec.split import clean_description, write_split_fasta


@pytest.fixture
def test_data_dir():
    return os.path.join(os.path.dirname(__file__), "test_data")


@pytest.fixture
def fasta_records(test_data_dir):
    """Fixture to parse records from multiple FASTA files into a dictionary."""
    files = [
        "multiple_records.fasta",
        "single_record.fasta",
        "has_empty_record.fasta",
        "has_empty_description.fasta",
    ]
    record_dict = {}
    for filename in files:
        file_path = os.path.join(test_data_dir, filename)
        with open(file_path, "r", encoding="utf-8") as input_file:
            records = list(SeqIO.parse(input_file, "fasta"))
            record_dict[filename] = records
    return record_dict


@pytest.mark.parametrize(
    "description, expected",
    [
        (
            'BBa_K620001_P_22737_Coding_"WT-F87A_(p450)"',
            "BBa_K620001_P_22737_Coding_WT-F87A_p450",
        ),
        ("long description" * 20, "longdescription" * 10),
        ("", ""),
    ],
)
def test_clean_description(description, expected):
    assert clean_description(description) == expected


@pytest.mark.parametrize(
    "filename",
    [
        "multiple_records.fasta",
        "single_record.fasta",
        "has_empty_record.fasta",
        "has_empty_description.fasta",
    ],
)
@patch("builtins.open", new_callable=mock_open)
@patch("os.path.join", side_effect=lambda a, b: f"{a}/{b}")
@patch("commec.split.SeqIO.parse")
def test_write_split_fasta(
    mock_seqio_parse,
    mock_os_path_join,
    mock_open,
    filename,
    test_data_dir,
    fasta_records,
):
    filepath = os.path.join(test_data_dir, filename)
    records = fasta_records[filename]
    mock_seqio_parse.return_value = records
    write_split_fasta(filepath)

    # Check the correct number of output files were opened (one input + as many outputs as records)
    assert mock_open.call_count == len(records) + 1

    for record in records:
        desc = clean_description(record.description)

        if desc:
            output_filename = f"{desc}.fasta"
        else:
            output_filename = f"{os.path.splitext(filename)[0]}-split-0.fasta"

        mock_os_path_join.assert_any_call(os.path.dirname(filepath), output_filename)
        mock_open.assert_any_call(
            os.path.join(os.path.dirname(filepath), output_filename),
            "w",
            encoding="utf-8",
        )
        mock_open().write.assert_any_call(f">{desc}{os.linesep}")
        mock_open().write.assert_any_call(f"{record.seq}")
