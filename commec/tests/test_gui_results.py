"""Result aggregation for the Commec GUI."""

import json

import pytest

from commec.gui import server


def _summary(tmp_path, statuses):
    output = {
        "queries": {
            f"query_{index}": {
                "length": 100,
                "status": {"screen_status": status},
            }
            for index, status in enumerate(statuses)
        }
    }
    (tmp_path / "test.output.json").write_text(json.dumps(output), encoding="utf-8")
    return server._summarize_output(tmp_path)


def test_all_too_short_run_retains_too_short_verdict(tmp_path):
    assert _summary(tmp_path, ["Skip (too short)"])["overall"] == "Skip (too short)"


@pytest.mark.parametrize(
    ("other_status", "expected"),
    [
        ("Pass", "Pass"),
        ("Warning", "Warning"),
        ("Skip (too long)", "Skip (too long)"),
        ("Skip", "Skip"),
        ("Flag", "Flag"),
    ],
)
def test_too_short_does_not_mask_actionable_status(tmp_path, other_status, expected):
    summary = _summary(tmp_path, ["Skip (too short)", other_status])
    assert summary["overall"] == expected
