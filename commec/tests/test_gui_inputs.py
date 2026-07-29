"""Input handling for the Commec GUI."""

import pytest
import yaml

from commec.gui import server


@pytest.fixture
def client(tmp_path, monkeypatch):
    monkeypatch.setitem(server.CFG, "password_hash", None)
    monkeypatch.setitem(server.CFG, "require_local_auth", False)
    monkeypatch.setitem(server.CFG, "runs_dir", tmp_path)
    monkeypatch.setitem(server.CFG, "default_databases", "")
    monkeypatch.setitem(
        server.CFG,
        "presets",
        [{"id": "test", "label": "Test", "config": {}}],
    )
    monkeypatch.setattr(server, "_missing_databases", lambda _preset: [])
    monkeypatch.setattr(server.threading.Thread, "start", lambda _thread: None)
    server.JOBS.clear()
    server.app.config.update(TESTING=True)
    yield server.app.test_client()
    server.JOBS.clear()


def _submit(client, **overrides):
    data = {
        "preset": "test",
        "label": "input-test",
        "sequence_text": ">short\nACGT\n",
    }
    data.update(overrides)
    return client.post(
        "/screen",
        data=data,
        environ_overrides={"REMOTE_ADDR": "127.0.0.1"},
    )


def test_short_sequences_are_passed_to_commec(tmp_path, client):
    """The GUI leaves sequence-length handling to the screening pipeline."""
    response = _submit(
        client,
        # Older clients may still submit the removed field. It must not
        # restore the GUI's former length filter.
        skip_short="1",
    )

    assert response.status_code == 200
    run_dir = next(path for path in tmp_path.iterdir() if path.is_dir())
    assert (run_dir / "input.fasta").read_text(encoding="utf-8") == ">short\nACGT\n"
    config = yaml.safe_load((run_dir / "config.used.yaml").read_text(encoding="utf-8"))
    assert config["databases"]["control_lists"]["regions"] == "all"


def test_selected_regions_are_recorded_in_run_config(tmp_path, client):
    response = _submit(client, regions="US,CA")

    assert response.status_code == 200
    run_dir = next(path for path in tmp_path.iterdir() if path.is_dir())
    config = yaml.safe_load((run_dir / "config.used.yaml").read_text(encoding="utf-8"))
    assert config["databases"]["control_lists"]["regions"] == "US,CA"


def test_unknown_region_is_rejected(tmp_path, client):
    response = _submit(client, regions="NOT_A_REGION")

    assert response.status_code == 400
    assert response.get_json() == {
        "error": "Unknown regulatory jurisdiction: NOT_A_REGION."
    }
    assert list(tmp_path.iterdir()) == []


def test_config_exposes_iso_region_choices(client):
    response = client.get(
        "/config",
        environ_overrides={"REMOTE_ADDR": "127.0.0.1"},
    )

    assert response.status_code == 200
    regions = response.get_json()["regions"]
    assert {"code": "US", "name": "United States", "group": False} in regions
