"""Input handling for the Commec GUI."""

import json

import pytest
import yaml

from commec.gui import server


@pytest.fixture
def client(tmp_path, monkeypatch):
    runs_dir = tmp_path / "runs"
    runs_dir.mkdir()
    databases = tmp_path / "databases"
    control_lists = databases / "control_lists"
    control_lists.mkdir(parents=True)
    (control_lists / "region_definitions.json").write_text(
        json.dumps(
            [
                {
                    "name": "Australia Group",
                    "acronym": "AG",
                    "regions": ["AU", "CA", "US"],
                },
                {
                    "name": "United Kingdom",
                    "acronym": "UK",
                    "regions": ["GB"],
                },
            ]
        ),
        encoding="utf-8",
    )
    direct_list = control_lists / "direct"
    direct_list.mkdir()
    (direct_list / "list_info.csv").write_text(
        "region_code\nBR\n",
        encoding="utf-8",
    )

    monkeypatch.setitem(server.CFG, "password_hash", None)
    monkeypatch.setitem(server.CFG, "require_local_auth", False)
    monkeypatch.setitem(server.CFG, "runs_dir", runs_dir)
    monkeypatch.setitem(server.CFG, "default_databases", str(databases))
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
    run_dir = next(path for path in (tmp_path / "runs").iterdir() if path.is_dir())
    assert (run_dir / "input.fasta").read_text(encoding="utf-8") == ">short\nACGT\n"
    config = yaml.safe_load((run_dir / "config.used.yaml").read_text(encoding="utf-8"))
    assert config["databases"]["control_lists"]["regions"] == "all"


def test_selected_regions_are_recorded_in_run_config(tmp_path, client):
    response = _submit(client, regions="US,CA")

    assert response.status_code == 200
    run_dir = next(path for path in (tmp_path / "runs").iterdir() if path.is_dir())
    config = yaml.safe_load((run_dir / "config.used.yaml").read_text(encoding="utf-8"))
    assert config["databases"]["control_lists"]["regions"] == "US,CA"


def test_hidden_single_country_alias_is_accepted(tmp_path, client):
    response = _submit(client, regions="UK")

    assert response.status_code == 200
    run_dir = next(path for path in (tmp_path / "runs").iterdir() if path.is_dir())
    config = yaml.safe_load((run_dir / "config.used.yaml").read_text(encoding="utf-8"))
    assert config["databases"]["control_lists"]["regions"] == "UK"


@pytest.mark.parametrize("region", ["NOT_A_REGION", "AQ"])
def test_unsupported_region_is_rejected(tmp_path, client, region):
    response = _submit(client, regions=region)

    assert response.status_code == 400
    assert response.get_json() == {
        "error": f"Unknown regulatory jurisdiction: {region}."
    }
    assert list((tmp_path / "runs").iterdir()) == []


def test_config_exposes_supported_region_choices(client):
    response = client.get(
        "/config",
        environ_overrides={"REMOTE_ADDR": "127.0.0.1"},
    )

    assert response.status_code == 200
    regions = response.get_json()["regions"]
    assert {
        "code": "US",
        "name": "United States",
        "group": False,
        "memberships": ["Australia Group"],
    } in regions
    assert {
        "code": "BR",
        "name": "Brazil",
        "group": False,
        "memberships": [],
    } in regions
    assert not any(region["code"] == "AQ" for region in regions)
    assert [region for region in regions if region["code"] == "AG"] == [
        {"code": "AG", "name": "Australia Group", "group": True}
    ]
    assert not any(region["code"] == "UK" for region in regions)
    assert {
        "code": "GB",
        "name": "United Kingdom",
        "group": False,
        "memberships": [],
    } in regions
