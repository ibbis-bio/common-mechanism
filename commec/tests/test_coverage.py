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
from commec.screeners.check_benign import _calculate_coverage
from commec.config.result import MatchRange

# Test the following hmmer configuration:
# 10-----------------50 (Largest, should stay.)
#               40-------60 (Extends 1. Should stay.)
#      20-------40      (encapsulated, but high score, should stay!)
# 10-----------------50 (lower score than 1, should be removed.)
#      20---30 (Different query, removed)
# 10------------------------------90 (different query)

# Example DataFrame
example_hmmer_01 = pd.DataFrame({
    "q. start": [100,  50,   0,   0],
    "q. end":   [200, 150, 100, 200],
})

# Example DataFrame
example_hmmer_01_output = pd.DataFrame({
    "q. start": [100,  50,   0,   0],
    "q. end":   [200, 150, 100, 200],
    "coverage_nt": [100, 50, 0,  100],
    "coverage_ratio": [1.0, 0.5, 0.0, 1.0]
})

reg_range_01 = MatchRange(0.0, 100, 200, 100, 200)

@pytest.mark.parametrize(
    "input_hmmer, input_region, expected_output_hmmer",
    [
        (example_hmmer_01, reg_range_01, example_hmmer_01_output),
    ]
)
def test_coverage_overlaps(
    input_hmmer : pd.DataFrame,
    input_region : MatchRange,
    expected_output_hmmer : pd.DataFrame
):
    """
    Checks common configurations that require trimming in Hmmer outputs,
    In particular partial overlaps, full encapsulations, score differences, and different queries.
    """
    trimmed_input = _calculate_coverage(input_hmmer, input_region)
    print("INPUT:")
    print(input_hmmer)
    print("TRIMMED:")
    print(trimmed_input)
    print("CORRECT:")
    print(expected_output_hmmer)
    assert trimmed_input.equals(expected_output_hmmer)
