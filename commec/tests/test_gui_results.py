"""Result aggregation for the Commec GUI."""

import json

import pytest

from commec.config.result import ScreenStatus
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
    summary = _summary(tmp_path, ["Skip (too short)"])
    assert summary["overall"] == "Skip (too short)"
    assert summary["overall_score"] == ScreenStatus.SKIP_SHORT.importance * 10
    assert summary["score_max"] == max(status.importance for status in ScreenStatus) * 10
    assert summary["queries"][0]["score"] == summary["overall_score"]


@pytest.mark.parametrize(
    ("statuses", "expected"),
    [
        (["Skip (too short)", "Pass"], "Pass"),
        (["Skip (too short)", "Warning"], "Warning"),
        (["Skip (too short)", "Skip (too long)"], "Skip (too long)"),
        (["Skip (too short)", "Skip"], "Skip (too short)"),
        (["Skip (too short)", "Flag"], "Flag"),
        (["Flag", "Incomplete"], "Incomplete"),
        (["Incomplete", "Error"], "Error"),
        (["Pass", "Future status"], "Future status"),
    ],
)
def test_run_uses_package_status_importance(tmp_path, statuses, expected):
    summary = _summary(tmp_path, statuses)
    assert summary["overall"] == expected


def test_recent_runs_refresh_cached_scores_from_package(tmp_path, monkeypatch):
    runs_dir = tmp_path / "runs"
    run_dir = runs_dir / "cached-run"
    run_dir.mkdir(parents=True)
    output = {
        "queries": {
            "query_1": {
                "length": 100,
                "status": {"screen_status": "Error"},
            }
        }
    }
    (run_dir / "test.output.json").write_text(json.dumps(output), encoding="utf-8")
    stale_meta = {
        "id": "cached-run",
        "label": "Cached run",
        "status": "done",
        "finished": 1,
        "summary": {"n": 1, "overall": "Pass"},
    }
    (run_dir / "meta.json").write_text(json.dumps(stale_meta), encoding="utf-8")

    monkeypatch.setitem(server.CFG, "runs_dir", runs_dir)
    monkeypatch.setitem(server.CFG, "password_hash", None)
    monkeypatch.setitem(server.CFG, "require_local_auth", False)
    server.JOBS.clear()
    response = server.app.test_client().get("/runs")

    assert response.status_code == 200
    run = response.get_json()["runs"][0]
    assert run["overall"] == "Error"
    assert run["overall_score"] == ScreenStatus.ERROR.importance * 10
    refreshed_meta = json.loads((run_dir / "meta.json").read_text(encoding="utf-8"))
    assert refreshed_meta["summary"]["overall"] == "Error"
