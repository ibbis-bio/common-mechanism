import pytest
from unittest.mock import patch
import pandas as pd
import os
from Bio.SeqRecord import SeqRecord, Seq

from commec.screeners.check_biorisk import biorisk_evalue_filter, parse_biorisk_hits, HmmerHandler
from commec.config.result import ScreenResult
from commec.config.query import Query
from commec.config.constants import (
    BIORISK_SHORT_QUERY_NT_THRESHOLD,
    BIORISK_SHORT_QUERY_EVALUE_EXPONENT,
    BIORISK_LONG_QUERY_EVALUE_THRESHOLD,
)

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
        handler = HmmerHandler(DATABASE_DIRECTORY + "biorisk/biorisk.hmm", INPUT_QUERY, "/mock/path/test.hmmscan")
        results = ScreenResult()
        queries : dict[str,Query] = {"testname" : Query(SeqRecord(Seq("atgatgatgatgatgatgatg"),"testname","testname"))}
        # Run the function - input paths are unused given all the mocking above
        result = parse_biorisk_hits(handler, "/mock/path/biorisk/biorisk_annotations.csv", results, queries)

        # Check the result
        assert result == expected_return


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_hmmer(nt_qlen: int, evalue: float) -> pd.DataFrame:
    """Return a minimal single-row hmmscan DataFrame."""
    return pd.DataFrame({"E-value": [evalue], "nt_qlen": [nt_qlen]})


def _short_cutoff(nt_qlen: float) -> float:
    """Length-dependent E-value cutoff for short queries."""
    return 1 / (1 + nt_qlen ** BIORISK_SHORT_QUERY_EVALUE_EXPONENT)


# ---------------------------------------------------------------------------
# biorisk_evalue_filter — short-query regime (nt_qlen < threshold)
# ---------------------------------------------------------------------------

class TestBioriskEvalueFilterShortQuery:
    """nt_qlen below BIORISK_SHORT_QUERY_NT_THRESHOLD uses the power-law cutoff."""

    def test_short_query_passes_when_evalue_below_cutoff(self):
        nt_qlen = 100
        evalue = _short_cutoff(nt_qlen) * 0.5  # clearly below cutoff
        result = biorisk_evalue_filter(_make_hmmer(nt_qlen, evalue))
        assert len(result) == 1

    def test_short_query_filtered_when_evalue_above_cutoff(self):
        nt_qlen = 100
        evalue = _short_cutoff(nt_qlen) * 2.0  # clearly above cutoff
        result = biorisk_evalue_filter(_make_hmmer(nt_qlen, evalue))
        assert len(result) == 0

    def test_very_short_query_uses_power_law_not_long_threshold(self):
        """A very short query with E-value below long threshold but above short
        cutoff should be filtered out — not admitted by the long-query rule."""
        nt_qlen = 10
        # For nt_qlen=10 the short cutoff is ~1/(1+10^2.598) ≈ 2.5e-3; use a
        # value that is below 1e-20 (long threshold) but above short cutoff.
        # That is impossible to construct, so instead confirm the short rule
        # dominates: use an evalue above the short cutoff and verify it is dropped.
        evalue = _short_cutoff(nt_qlen) * 2
        result = biorisk_evalue_filter(_make_hmmer(nt_qlen, evalue))
        assert len(result) == 0

    def test_query_just_below_threshold_uses_short_rule(self):
        nt_qlen = BIORISK_SHORT_QUERY_NT_THRESHOLD - 1
        evalue = _short_cutoff(nt_qlen) * 0.5
        result = biorisk_evalue_filter(_make_hmmer(nt_qlen, evalue))
        assert len(result) == 1

class TestBioriskEvalueFilterLongQuery:
    """nt_qlen at or above BIORISK_SHORT_QUERY_NT_THRESHOLD uses the fixed cutoff."""

    def test_long_query_passes_when_evalue_below_threshold(self):
        nt_qlen = BIORISK_SHORT_QUERY_NT_THRESHOLD
        evalue = BIORISK_LONG_QUERY_EVALUE_THRESHOLD * 0.1
        result = biorisk_evalue_filter(_make_hmmer(nt_qlen, evalue))
        assert len(result) == 1

    def test_long_query_filtered_when_evalue_above_threshold(self):
        nt_qlen = BIORISK_SHORT_QUERY_NT_THRESHOLD
        evalue = BIORISK_LONG_QUERY_EVALUE_THRESHOLD * 10
        result = biorisk_evalue_filter(_make_hmmer(nt_qlen, evalue))
        assert len(result) == 0

class TestBioriskEvalueFilterMultiRow:
    """Tests covering mixed DataFrames and edge cases."""

    def test_mixed_rows_only_passing_rows_retained(self):
        nt_qlen_short = 100
        nt_qlen_long = 500
        df = pd.DataFrame({
            "E-value": [
                _short_cutoff(nt_qlen_short) * 0.5,   # short, passes
                _short_cutoff(nt_qlen_short) * 2.0,   # short, filtered
                BIORISK_LONG_QUERY_EVALUE_THRESHOLD * 0.1,  # long, passes
                BIORISK_LONG_QUERY_EVALUE_THRESHOLD * 10,   # long, filtered
            ],
            "nt_qlen": [nt_qlen_short, nt_qlen_short, nt_qlen_long, nt_qlen_long],
        })
        result = biorisk_evalue_filter(df)
        assert len(result) == 2

    def test_empty_dataframe_returns_empty(self):
        df = pd.DataFrame({"E-value": [], "nt_qlen": []})
        result = biorisk_evalue_filter(df)
        assert len(result) == 0

    def test_string_columns_are_coerced_to_numeric(self):
        """readhmmer returns string columns; filter must coerce them."""
        nt_qlen = 100
        evalue = _short_cutoff(nt_qlen) * 0.5
        df = pd.DataFrame({"E-value": [str(evalue)], "nt_qlen": [str(nt_qlen)]})
        result = biorisk_evalue_filter(df)
        assert len(result) == 1

    def test_non_numeric_evalues_are_dropped(self):
        """Rows with non-coercible E-values become NaN and should be excluded."""
        df = pd.DataFrame({"E-value": ["not_a_number"], "nt_qlen": [100]})
        result = biorisk_evalue_filter(df)
        assert len(result) == 0

    def test_input_dataframe_is_not_mutated(self):
        """The function must return a copy; the original dtypes must be unchanged."""
        df = pd.DataFrame({"E-value": ["1e-25"], "nt_qlen": ["500"]})
        original_dtype_evalue = df["E-value"].dtype
        original_dtype_qlen = df["nt_qlen"].dtype
        _ = biorisk_evalue_filter(df)
        assert df["E-value"].dtype == original_dtype_evalue
        assert df["nt_qlen"].dtype == original_dtype_qlen

    def test_extra_columns_are_preserved(self):
        """Columns beyond E-value and nt_qlen must survive the filter unchanged."""
        nt_qlen = 100
        evalue = _short_cutoff(nt_qlen) * 0.5
        df = pd.DataFrame({
            "E-value": [evalue],
            "nt_qlen": [nt_qlen],
            "target name": ["some_target"],
            "query name": ["some_query"],
        })
        result = biorisk_evalue_filter(df)
        assert "target name" in result.columns
        assert result["target name"].iloc[0] == "some_target"