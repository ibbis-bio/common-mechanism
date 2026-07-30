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


@pytest.mark.parametrize(
    ("run_status", "result_status"),
    [
        ("error", None),
        ("done", "Error"),
    ],
)
def test_errored_run_can_be_rerun_in_a_clean_directory(
    tmp_path, monkeypatch, run_status, result_status
):
    runs_dir = tmp_path / "runs"
    source = runs_dir / "failed-run"
    broken = source / "output_failed" / "broken.partial"
    broken.parent.mkdir(parents=True)
    broken.write_text("partial output", encoding="utf-8")
    (source / "input.fasta").write_text(">sample\nACGT\n", encoding="utf-8")
    (source / "config.used.yaml").write_text("threads: 1\n", encoding="utf-8")
    meta = {
        "id": source.name,
        "label": "Failed run",
        "prefix": "failed",
        "status": run_status,
        "created": 1,
        "finished": 2,
    }
    (source / "meta.json").write_text(json.dumps(meta), encoding="utf-8")
    if result_status:
        output = {
            "queries": {
                "sample": {
                    "length": 4,
                    "status": {"screen_status": result_status},
                }
            }
        }
        (source / "failed.output.json").write_text(
            json.dumps(output), encoding="utf-8"
        )

    monkeypatch.setitem(server.CFG, "runs_dir", runs_dir)
    monkeypatch.setitem(server.CFG, "password_hash", None)
    monkeypatch.setitem(server.CFG, "require_local_auth", False)
    monkeypatch.setitem(server.CFG, "default_databases", "")
    monkeypatch.setitem(server.CFG, "threads", None)
    monkeypatch.setattr(server.threading.Thread, "start", lambda _thread: None)
    server.JOBS.clear()
    client = server.app.test_client()

    detail = client.get(f"/results/{source.name}")
    assert detail.status_code == 200
    assert detail.get_json()["rerunnable"] is True
    recent = client.get("/runs").get_json()["runs"]
    assert recent[0]["rerunnable"] is True

    response = client.post(f"/runs/{source.name}/rerun")

    assert response.status_code == 200
    new_id = response.get_json()["job_id"]
    assert new_id != source.name
    job = server.JOBS[new_id]
    assert job["label"] == "Failed run (rerun)"
    assert job["rerun_of"] == source.name
    assert "-R" not in job["cmd"]
    assert str(source) not in job["cmd"]
    copied_files = {
        path.relative_to(job["dir"]).as_posix()
        for path in job["dir"].rglob("*")
        if path.is_file()
    }
    assert copied_files == {"config.used.yaml", "input.fasta"}
    assert broken.read_text(encoding="utf-8") == "partial output"
    server.JOBS.clear()


def test_successful_run_cannot_be_rerun(tmp_path, monkeypatch):
    runs_dir = tmp_path / "runs"
    source = runs_dir / "successful-run"
    source.mkdir(parents=True)
    (source / "input.fasta").write_text(">sample\nACGT\n", encoding="utf-8")
    (source / "config.used.yaml").write_text("threads: 1\n", encoding="utf-8")
    meta = {
        "id": source.name,
        "label": "Successful run",
        "prefix": "successful",
        "status": "done",
    }
    (source / "meta.json").write_text(json.dumps(meta), encoding="utf-8")
    _summary(source, ["Pass"])

    monkeypatch.setitem(server.CFG, "runs_dir", runs_dir)
    monkeypatch.setitem(server.CFG, "password_hash", None)
    monkeypatch.setitem(server.CFG, "require_local_auth", False)
    server.JOBS.clear()
    response = server.app.test_client().post(f"/runs/{source.name}/rerun")

    assert response.status_code == 400
    assert response.get_json()["error"] == "That run did not fail; nothing to rerun."
