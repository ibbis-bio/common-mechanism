"""Authentication boundaries for the Commec GUI."""

import pytest

from commec.gui import server


@pytest.fixture
def client(monkeypatch):
    monkeypatch.setitem(server.CFG, "password_hash", None)
    monkeypatch.setitem(server.CFG, "require_local_auth", False)
    server.app.config.update(TESTING=True)
    return server.app.test_client()


def test_hashless_gui_allows_localhost(client):
    response = client.get("/config", environ_overrides={"REMOTE_ADDR": "127.0.0.1"})

    assert response.status_code == 200


def test_hashless_gui_rejects_remote_access(client):
    response = client.get("/config", environ_overrides={"REMOTE_ADDR": "192.168.1.10"})

    assert response.status_code == 503
    assert response.get_json() == {
        "error": "Remote access is unavailable until an operator password is set."
    }
