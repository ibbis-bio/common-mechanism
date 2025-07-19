import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from commec.config.query import Query, QueryValueError

@pytest.fixture
def seq_record():
    """
    Fixture to generate a SeqRecord with 
    defined coding ('t') and non-coding ('a') regions.
    Total size = 150
    non-coding regions = 100
    """
    non_coding_1 = "a" * 50
    coding_1 = "t" * 20
    non_coding_2 = "a" * 30
    coding_2 = "t" * 30
    non_coding_3 = "a" * 20  # Final non-coding segment

    sequence = non_coding_1 + coding_1 + non_coding_2 + coding_2 + non_coding_3
    return SeqRecord(Seq(sequence), id="test_seq", description="")

@pytest.fixture
def non_coding_regions():
    """
    Fixture to generate non-coding region tuples based on the sequence definition.
    Uses the same lengths as in seq_record() to compute (start, end) values.
    """
    regions = [
        (1, 50),       # First 'a' region
        (71, 100),     # Second 'a' region (starts after first coding region)
        (131, 150)     # Third 'a' region (starts after second coding region)
    ]
    return regions

# 0 based coordinates:
# NT COORDS: 0 - 49    50 - 69   70 - 99   100 - 129   130 - 149
# NC COORDS: 0 - 49              50 - 79                80 -  99

# 1 based coordinates:
# NT COORDS: 1 - 50    51 - 70   71 - 100   101 - 130   131 - 150
# NC COORDS: 1 - 50              51 -  80                81 - 100

@pytest.fixture
def test_cases():
    """
    Fixture providing a list of (input_coordinate, expected_output) tuples.
    The input is a sequence coordinate, and the expected output is its transformed coordinate.
    """
    return [
        (1, 1),
        (10, 10),
        (50, 50),
        (51, 71),
        (80, 100),
        (81, 131),
        (99, 149),
        (100, 150),
    ]

def test_coordinate_conversion(seq_record, non_coding_regions, test_cases):
    """
    Placeholder test function for coordinate conversion.
    """
    # Query setup
    test_query : Query = Query(seq_record)
    test_query.non_coding_regions = non_coding_regions

    # Test Correct coords:
    for nc, nt in test_cases:
        assert nt == test_query.nc_to_nt_query_coords(nc)

    # Test Failure out of bounds.
    try:
        _x = test_query.nc_to_nt_query_coords(test_cases[-1][0]+1)
        assert False
    except QueryValueError:
        assert True

    # Test Failure out of bounds.
    try:
        _x = test_query.nc_to_nt_query_coords(0)
        assert False
    except QueryValueError:
        assert True

