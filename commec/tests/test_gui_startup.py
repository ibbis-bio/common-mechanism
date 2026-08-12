"""Startup diagnostics for the Commec GUI."""

import pytest

from commec.gui import server


@pytest.mark.parametrize(
    "failure",
    [SystemExit(1), PermissionError("Permission denied")],
)
def test_privileged_port_failure_prints_hint(monkeypatch, capsys, failure):
    def fail_to_start(**_kwargs):
        raise failure

    monkeypatch.setattr(server.app, "run", fail_to_start)

    with pytest.raises(type(failure)):
        server._serve("127.0.0.1", 443, None)

    error = capsys.readouterr().err
    assert "Port 443 is privileged on many systems" in error
    assert "commec gui --port 8765" in error


def test_high_port_failure_does_not_print_privileged_port_hint(monkeypatch, capsys):
    def fail_to_start(**_kwargs):
        raise SystemExit(1)

    monkeypatch.setattr(server.app, "run", fail_to_start)

    with pytest.raises(SystemExit):
        server._serve("127.0.0.1", 8765, None)

    assert "privileged" not in capsys.readouterr().err
