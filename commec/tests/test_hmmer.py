"""
Simple tests to ensure the correct trimming of overlapping components in a Hmmer database
when parsed to remove_overlaps.
The following behaviour is expected:
 * Fully encapsulated hits, should be removed.
 * Partially overlapping hits, should both be kept (To maximise extents)
 * Hits from different queries are independant in logic.
"""

import pandas as pd
import pytest
from commec.tools.hmmer import remove_overlaps, get_frame_from_query_name

# Test the following hmmer configuration:
# 10-----------------50 (Largest, should stay.)
#               40-------60 (Extends 1. Should stay.)
#      20-------40      (encapsulated, but high score, should stay!)
# 10-----------------50 (lower score than 1, should be removed.)
#      20---30 (Different query, removed)
# 10------------------------------90 (different query)

# Example DataFrame
example_hmmer_01 = pd.DataFrame({
    "query name": ["one","one","one","one","two", "two"],
    "q. start": [10, 40, 20, 10, 20, 10],
    "q. end":   [50, 60, 40, 50, 30, 90],
    "score":    [3, 5, 6, 1, 1, 2]
})

# Example DataFrame
example_hmmer_01_output = pd.DataFrame({
    "query name": ["one","one", "one", "two"],
    "q. start": [10, 40, 20, 10],
    "q. end":   [50, 60, 40, 90],
    "score":    [3, 5, 6, 2]
})

@pytest.mark.parametrize(
    "input_hmmer, expected_output_hmmer",
    [
        (example_hmmer_01, example_hmmer_01_output),
    ]
)
def test_remove_overlaps(
    input_hmmer : pd.DataFrame,
    expected_output_hmmer : pd.DataFrame
):
    """
    Checks common configurations that require trimming in Hmmer outputs,
    In particular partial overlaps, full encapsulations, score differences, and different queries.
    """
    trimmed_input = remove_overlaps(input_hmmer)
    print("INPUT:")
    print(input_hmmer)
    print("TRIMMED:")
    print(trimmed_input)
    print("CORRECT:")
    print(expected_output_hmmer)
    assert trimmed_input.equals(expected_output_hmmer)


def test_get_frame_from_query_name():
    example_hmmer_01_output = pd.DataFrame({
        "query name": ["one","one", "one", "two"],
        "q. start": [10, 40, 20, 10],
        "q. end":   [50, 60, 40, 90],
        "score":    [3, 5, 6, 2]
    })
    get_frame_from_query_name(example_hmmer_01_output["query name"])