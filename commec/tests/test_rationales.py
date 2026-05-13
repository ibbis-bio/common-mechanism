"""
Unit tests for controlled rationale outcomes
"""

from commec.tests.screen_factory import (
    ScreenTesterFactory,
    ScreenStep
)
from commec.config.result import ScreenStatus, Rationale
from commec.config.constants import MINIMUM_QUERY_LENGTH, MAXIMUM_QUERY_LENGTH


def test_hmmer(tmp_path):
    """
    When there are hits to Biorisk with a large E-value, but no other hits, and we 
     are running in the skip taxonomy mode, we correctly label the outcome
     as warning, however the rationale is set to "Matches to ." instead of
     the correct Rationale text indicating no hits.
    """
    screen_test = ScreenTesterFactory("low_evalue_hmmer", tmp_path)
    screen_test.add_query("query1",1200)
    screen_test.add_hit(ScreenStep.BIORISK, "query1", 100, 200, "HighEvalueHit", "HEH", 500, regulated=True, evalue = 100.0)
    result = screen_test.run("--skip-tx")
    assert result.queries["query1"].status.screen_status == ScreenStatus.WARN
    assert result.queries["query1"].status.rationale == str(Rationale.NO_HITS_SKIP_NOTE)

def test_skip_rationales(tmp_path):
    """
    Ensure when a skip occurs due to query length, 
    the rationale and screen status are correct. 
    """
    screen_test = ScreenTesterFactory("low_evalue_hmmer", tmp_path)
    screen_test.add_query("query1",MINIMUM_QUERY_LENGTH - 1)
    screen_test.add_query("query2",MAXIMUM_QUERY_LENGTH + 1)
    result = screen_test.run()
    assert result.queries["query1"].status.screen_status == ScreenStatus.SKIP_SHORT
    assert result.queries["query2"].status.screen_status == ScreenStatus.SKIP_LONG
    assert result.queries["query1"].status.rationale == str(Rationale.SKIPPED) + " " + str(Rationale.TOO_SHORT)
    assert result.queries["query2"].status.rationale == str(Rationale.SKIPPED) + " " + str(Rationale.TOO_LONG)
