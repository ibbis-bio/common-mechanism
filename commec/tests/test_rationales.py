"""
Unit tests for controlled rationale outcomes
"""

from commec.tests.screen_factory import (
    ScreenTesterFactory,
    ScreenStep
)
from commec.config.result import ScreenStatus, Rationale

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
    assert result.queries["query1"].status.rationale == str(Rationale.NOTHING)
