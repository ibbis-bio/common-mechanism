"""
Simple tests to ensure the correct trimming of overlapping components in a Hmmer database
when parsed to remove_overlaps.
The following behaviour is expected:
 * Fully encapsulated hits, should be removed.
 * Partially overlapping hits, should both be kept (To maximise extents)
 * Hits from different queries are independant in logic.
"""

from unittest.mock import patch

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal, assert_series_equal
from commec.tools.hmmer import (
    remove_overlaps,
    get_frame_from_query_name,
    set_query_nt_coordinates,
)

# First test - multiple queries and overlaps
# Test the following hmmer configuration:
# 10-----------------50 (Largest, should stay.)
#               40-------60 (Extends 1. Should stay.)
#      20-------40      (encapsulated, but high score, should stay!)
# 10-----------------50 (lower score than 1, should be removed.)
#      20---30 (Different query, removed)
# 10------------------------------90 (different query)
hmmer_with_overlaps = pd.DataFrame(
    {
        "query name": ["one", "one", "one", "one", "two", "two"],
        "q. start": [10, 40, 20, 10, 20, 10],
        "q. end": [50, 60, 40, 50, 30, 90],
        "score": [3, 5, 6, 1, 1, 2],
    }
)
expected_hmmer = pd.DataFrame(
    {
        "query name": ["one", "one", "one", "two"],
        "q. start": [10, 40, 20, 10],
        "q. end": [50, 60, 40, 90],
        "score": [3, 5, 6, 2],
    }
)

# Second test - two hits with equal scores
hmmer_with_equal_overlaps = pd.DataFrame(
    {
        "query name": ["one", "one"],
        "q. start": [10, 10],
        "q. end": [50, 50],
        "score": [5, 5],
    }
)
expected_hmmer_with_equal = pd.DataFrame(
    {"query name": ["one"], "q. start": [10], "q. end": [50], "score": [5]}
)

# Third test - only a single hit
hmmer_with_single_hit = pd.DataFrame(
    {"query name": ["one"], "q. start": [10], "q. end": [50], "score": [5]}
)


@pytest.mark.parametrize(
    "input_hmmer, expected_output_hmmer",
    [
        (hmmer_with_overlaps, expected_hmmer),
        (hmmer_with_equal_overlaps, expected_hmmer_with_equal),
        (hmmer_with_single_hit, hmmer_with_single_hit.copy()),
    ],
)
def test_remove_overlaps(
    input_hmmer: pd.DataFrame, expected_output_hmmer: pd.DataFrame
):
    """
    Checks common configurations that require trimming in Hmmer outputs,
    In particular partial overlaps, full encapsulations, score differences, and different queries.
    """
    # No need to actually print warnings from this
    with patch("commec.tools.hmmer.set_query_nt_coordinates"):
        trimmed_input = remove_overlaps(input_hmmer)
        print("INPUT:")
        print(input_hmmer)
        print("TRIMMED:")
        print(trimmed_input)
        print("CORRECT:")
        print(expected_output_hmmer)
        assert trimmed_input.equals(expected_output_hmmer)


def test_get_frame_from_query_name():
    # Assumes frame 1 if it can't parse an int between 1 and 6 from the name
    cant_parse = pd.Series(["one", "1_F", "R_146810", "KL8700289", ""])
    with pytest.warns(UserWarning, match="Could not parse frame"):
        assert_series_equal(
            pd.Series([1, 1, 1, 1, 1]), get_frame_from_query_name(cant_parse)
        )

    can_parse = pd.Series(["F1", "F_2", "R_1468106", "1_F_5"])
    assert_series_equal(pd.Series([1, 2, 6, 5]), get_frame_from_query_name(can_parse))


def test_set_query_coordinates():
    """
    Diagram showing six-frame translation of a 10 amino acid / 30 nt sequence:

        1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
        [   1  ] [   2  ] [   3  ] [   4  ] [   5  ] [   6  ] [   7  ] [   8  ] [   9  ] [  10  ]
           [  1  ]  [  2  ]  [  3  ] [   4  ] [   5  ] [   6  ] [   7  ] [   8  ] [   9  ]
              [  1  ]  [  2  ]  [  3  ] [   4  ] [   5  ] [   6  ] [   7  ] [   8  ] [   9  ]
        [  10  ] [   9  ] [   8  ] [   7  ] [   6  ] [   5  ] [   4  ] [   3  ] [   2  ] [   1  ]
          [   9  ] [   8  ] [   7  ] [   6  ] [   5  ] [   4  ] [   3  ] [   2  ] [   1  ]
              [   9  ] [   8  ] [   7  ] [   6  ] [   5  ] [   4  ] [   3  ] [   2  ] [   1  ]
    """
    input_hmmer = pd.DataFrame(
        {
            "query name": ["F1", "F2", "F3", "R1", "R2", "R3", "R3"],
            "frame": [1, 2, 3, 4, 5, 6, 6],
            "ali from": [1, 2, 3, 1, 2, 3, 4],
            "ali to": [4, 5, 6, 4, 5, 6, 4],
            "qlen": [10, 9, 9, 10, 9, 9, 9],
            "nt len": [30, 30, 30, 30, 30, 30, 30],
        }
    )

    hmmer_with_query_coords = pd.DataFrame(
        {
            "query name": ["F1", "F2", "F3", "R1", "R2", "R3", "R3"],
            "frame": [1, 2, 3, 4, 5, 6, 6],
            "ali from": [1, 2, 3, 1, 2, 3, 4],
            "ali to": [4, 5, 6, 4, 5, 6, 4],
            "qlen": [10, 9, 9, 10, 9, 9, 9],
            "nt len": [30, 30, 30, 30, 30, 30, 30],
            "q. start": [1, 5, 9, 19, 14, 12, 18],
            "q. end": [12, 16, 20, 30, 25, 23, 20],
        }
    )

    set_query_nt_coordinates(input_hmmer)
    print("PROCESSED:")
    print(input_hmmer)
    print("CORRECT:")
    print(hmmer_with_query_coords)
    assert_frame_equal(input_hmmer, hmmer_with_query_coords)
