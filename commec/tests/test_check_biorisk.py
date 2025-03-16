import pytest
from unittest.mock import patch
import pandas as pd

from commec.screeners.check_biorisk import check_biorisk


@pytest.mark.parametrize(
    "annotations_exists, is_empty, has_hits, expected_return",
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
def test_check_biorisk_return_codes(annotations_exists, is_empty, has_hits, expected_return):
    mock_hit_df = pd.DataFrame(
        {
            "target name": ["test_id"],
            "E-value": [1e-30],
            "ali from": [100],
            "ali to": [200],
            "qlen": [1000],
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
        patch("commec.screeners.check_biorisk.trimhmmer", return_value=mock_hit_df),
        patch("commec.screeners.check_biorisk.HmmerHandler.is_empty", return_value=is_empty),
        patch("commec.screeners.check_biorisk.HmmerHandler.has_hits", return_value=has_hits),
    ):
        # Run the function - input paths are unused given all the mocking above
        result = check_biorisk("/mock/path/test.hmmscan", "/mock/path/biorisk_db")

        # Check the result
        assert result == expected_return
