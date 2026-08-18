"""Persisted age-based retention for Commec GUI runs."""

import json
import os

import pytest

from commec.gui import server


NOW = 2_000_000_000
DAY = 86400


@pytest.fixture
def client(tmp_path, monkeypatch):
    runs_dir = tmp_path / "runs"
    runs_dir.mkdir()
    monkeypatch.setitem(server.CFG, "runs_dir", runs_dir)
    monkeypatch.setitem(server.CFG, "retention_days", 0)
    monkeypatch.setitem(server.CFG, "password_hash", None)
    monkeypatch.setitem(server.CFG, "require_local_auth", False)
    server.JOBS.clear()
    server.app.config.update(TESTING=True)
    yield server.app.test_client()
    server.JOBS.clear()


def _run_dir(runs_dir, name, status, finished):
    path = runs_dir / name
    path.mkdir()
    (path / "meta.json").write_text(
        json.dumps({
            "id": name,
            "label": name,
            "status": status,
            "finished": finished,
        }),
        encoding="utf-8",
    )
    return path


def test_retention_defaults_to_forever(client):
    assert server.CFG["retention_days"] == 0
    assert server._load_retention_days() == 0
    assert client.get("/config").get_json()["retention_days"] == 0


def test_retention_setting_is_persisted_and_exposed(client):
    response = client.post("/retention", json={"days": 30})

    assert response.status_code == 200
    assert response.get_json() == {
        "ok": True,
        "retention_days": 30,
        "deleted": 0,
    }
    assert json.loads(server._retention_path().read_text(encoding="utf-8")) == {
        "days": 30,
    }
    assert server._load_retention_days() == 30
    assert client.get("/config").get_json()["retention_days"] == 30


def test_setting_retention_immediately_deletes_only_expired_runs(
    client, monkeypatch
):
    runs_dir = server.CFG["runs_dir"]
    expired = _run_dir(runs_dir, "expired", "done", NOW - 31 * DAY)
    retained = _run_dir(runs_dir, "retained", "error", NOW - 29 * DAY)
    boundary = _run_dir(runs_dir, "boundary", "interrupted", NOW - 30 * DAY)
    monkeypatch.setattr(server.time, "time", lambda: NOW)

    response = client.post("/retention", json={"days": 30})

    assert response.status_code == 200
    assert response.get_json()["deleted"] == 1
    assert not expired.exists()
    assert retained.exists()
    assert boundary.exists()


def test_sweep_protects_live_and_nonterminal_runs(client):
    runs_dir = server.CFG["runs_dir"]
    live = _run_dir(runs_dir, "live", "running", NOW - 10 * DAY)
    unknown = _run_dir(runs_dir, "unknown", "starting", NOW - 10 * DAY)
    server.CFG["retention_days"] = 1
    server.JOBS["live"] = {"dir": live, "done": False}

    deleted = server._sweep_runs(now=NOW)

    assert deleted == []
    assert live.exists()
    assert unknown.exists()


def test_finished_timestamp_takes_precedence_over_directory_mtime(client):
    runs_dir = server.CFG["runs_dir"]
    run = _run_dir(runs_dir, "old-finish", "done", NOW - 10 * DAY)
    os.utime(run, (NOW, NOW))
    server.CFG["retention_days"] = 1

    assert server._sweep_runs(now=NOW) == ["old-finish"]
    assert not run.exists()


@pytest.mark.parametrize("value", [-1, 1.5, True, None, "abc", "01"])
def test_invalid_retention_is_rejected(client, value):
    response = client.post("/retention", json={"days": value})

    assert response.status_code == 400
    assert response.get_json() == {
        "error": "Retention days must be a non-negative whole number."
    }
    assert server.CFG["retention_days"] == 0
    assert not server._retention_path().exists()
