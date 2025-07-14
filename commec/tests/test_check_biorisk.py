import pytest
from unittest.mock import patch
import pandas as pd
import os
from Bio.SeqRecord import SeqRecord, Seq

from commec.screeners.check_biorisk import update_biorisk_data_from_database, HmmerHandler
from commec.config.result import ScreenResult
from commec.config.query import Query

INPUT_QUERY = os.path.join(os.path.dirname(__file__), "test_data/single_record.fasta")
DATABASE_DIRECTORY = os.path.join(os.path.dirname(__file__), "test_dbs/")

@pytest.mark.parametrize(
    "annotations_exists, has_empty_output, has_hits, expected_return",
    [
        # Case 1: annotations file doesn't exist
        (False, False, False, 1),
        # Case 2: HMMER output is empty or doesn't exist
        (True, True, False, 1),
        # Case 3: No hits detected (successful pass)
        (True, False, False, 0),
        # Case 4: Successful execution with hits
        (True, False, True, 0),
    ],
)
def test_check_biorisk_return_codes(annotations_exists, has_empty_output, has_hits, expected_return):
    mock_hit_df = pd.DataFrame(
        {
            "target name": ["test_id"],
            "query name": ["testname_1"],
            "E-value": [1e-30],
            "ali from": [100],
            "ali to": [200],
            "qlen": [1000],
            "frame" : 1
        }
    )

    mock_annot_df = pd.DataFrame(
        {"ID": ["test_id"], "description": ["test description"], "Must flag": [True]}
    )

    # No filesystem interactions, patch ALL the things
    with (
        patch("os.path.exists", return_value=annotations_exists),
        patch("pandas.read_csv", return_value=mock_annot_df),
        patch("commec.screeners.check_biorisk.readhmmer", return_value=mock_hit_df),
        patch("commec.screeners.check_biorisk.remove_overlaps", return_value=mock_hit_df),
        patch("commec.screeners.check_biorisk.HmmerHandler.has_empty_output", return_value=has_empty_output),
        patch("commec.screeners.check_biorisk.HmmerHandler.has_hits", return_value=has_hits),
    ):
        handler = HmmerHandler(DATABASE_DIRECTORY + "biorisk_db/biorisk.hmm", INPUT_QUERY, "/mock/path/test.hmmscan")
        results = ScreenResult()
        queries : dict[str,Query] = {"testname" : Query(SeqRecord(Seq("atgatgatgatgatgatgatg"),"testname","testname"))}
        # Run the function - input paths are unused given all the mocking above
        result = update_biorisk_data_from_database(handler, "/mock/path/biorisk_db/biorisk_annotations.csv", results, queries)

        # Check the result
        assert result == expected_return