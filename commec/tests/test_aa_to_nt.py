"""
Simple tests to ensure the correct trimming of overlapping components in a Hmmer database
when parsed to remove_overlaps.
The following behaviour is expected:
 * Fully encapsulated hits, should be removed.
 * Partially overlapping hits, should both be kept (To maximise extents)
 * Hits from different queries are independant in logic.
"""

import pandas as pd
from pandas.testing import assert_frame_equal
import pytest
from commec.tools.hmmer import recalculate_hmmer_query_coordinates

# Example DataFrame
example_hmmer_01 = pd.DataFrame({
    "query name": ["F1","F2","F3","R1","R2", "R3"],
    "frame":    [1,2,3,4,5,6],
    "ali from": [1, 2, 3, 1, 2, 3],
    "ali to":   [4, 5, 6, 4, 5, 6],
    "nt_qlen":     [31, 31, 31, 31, 31, 31]
})

# Logically, we would expect this to match Fwd and Rev all frame AA to NT coordinates:
# 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
# [  1  ]  [  2  ]  [  3  ] [   4  ] [   5  ] [   6  ] [   7  ] [   8  ] [   9  ] [  10  ]
#    [  1  ]  [  2  ]  [  3  ] [   4  ] [   5  ] [   6  ] [   7  ] [   8  ] [   9  ] [  10  ]
#       [  1  ]  [  2  ]  [  3  ] [   4  ] [   5  ] [   6  ] [   7  ] [   8  ] [   9  ]
# [  10 ]  [  9  ]  [  8  ] [   7  ] [   6  ] [   5  ] [   4  ] [   3  ] [   2  ] [   1  ]
#       [  9  ]  [  8  ] [   7  ] [   6  ] [   5  ] [   4  ] [   3  ] [   2  ] [   1  ]
#    [  9  ]  [  8  ] [   7  ] [   6  ] [   5  ] [   4  ] [   3  ] [   2  ] [   1  ] [  -1  ]

# However, HMMER Biorisk behaves like the following 
#(According to testing with BBa_I766605 YopH-EE under medium constitutive promotor, and reverse complements.)
# 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
# [  1  ]  [  2  ]  [  3  ] [   4  ] [   5  ] [   6  ] [   7  ] [   8  ] [   9  ] [  10  ]
#    [  1  ]  [  2  ]  [  3  ] [   4  ] [   5  ] [   6  ] [   7  ] [   8  ] [   9  ] [  10  ]
#       [  1  ]  [  2  ]  [  3  ] [   4  ] [   5  ] [   6  ] [   7  ] [   8  ] [   9  ]
#       [  10 ]  [  9  ]  [  8  ] [   7  ] [   6  ] [   5  ] [   4  ] [   3  ] [   2  ] [   1  ]
#             [  9  ]  [  8  ] [   7  ] [   6  ] [   5  ] [   4  ] [   3  ] [   2  ] [   1  ]
#          [  9  ]  [  8  ] [   7  ] [   6  ] [   5  ] [   4  ] [   3  ] [   2  ] [   1  ] [  -1  ]

# Example DataFrame
original_expected_output = pd.DataFrame({
    "query name": ["F1","F2","F3","R1","R2", "R3"],
    "frame":    [1,2,3,4,5,6],
    "ali from": [1, 2, 3, 1, 2, 3],
    "ali to":   [4, 5, 6, 4, 5, 6],
    "nt_qlen":  [31, 31, 31, 31, 31, 31],
    "q. start": [ 1,  5,  9,  19,  15, 11],
    "q. end":   [12, 16, 20,  30,  26, 22]
})

# Example DataFrame which matches biorisk hmmer outputs:
example_hmmer_01_output = pd.DataFrame({
    "query name": ["F1","F2","F3","R1","R2", "R3"],
    "frame":    [1,2,3,4,5,6],
    "ali from": [1, 2, 3, 1, 2, 3],
    "ali to":   [4, 5, 6, 4, 5, 6],
    "nt_qlen":  [31, 31, 31, 31, 31, 31],
    "q. start": [ 1,  5,  9,  21,  17, 13],
    "q. end":   [12, 16, 20,  32,  28, 24]
})

@pytest.mark.parametrize(
    "input_hmmer, expected_output_hmmer",
    [
        (example_hmmer_01, example_hmmer_01_output),
    ]
)
def test_hmmer_overlaps(
    input_hmmer : pd.DataFrame,
    expected_output_hmmer : pd.DataFrame
):
    """
    Checks common configurations that require trimming in Hmmer outputs,
    In particular partial overlaps, 
    full encapsulations, score differences, 
    and different queries.
    """
    print("INPUT:")
    print(input_hmmer)
    print(input_hmmer.dtypes)
    recalculate_hmmer_query_coordinates(input_hmmer)
    print("PROCESSED:")
    print(input_hmmer)
    print(input_hmmer.dtypes)
    print("CORRECT:")
    print(expected_output_hmmer)
    print(expected_output_hmmer.dtypes)
    assert_frame_equal(input_hmmer, expected_output_hmmer)
    #assert input_hmmer.equals(expected_output_hmmer)
