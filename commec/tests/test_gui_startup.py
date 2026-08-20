"""Startup diagnostics for the Commec GUI."""

from pathlib import Path

import pytest

from commec.gui import server


def test_deployment_presets_override_missing_packaged_presets(monkeypatch, tmp_path):
    deployed = tmp_path / "presets.yaml"
    deployed.write_text("presets: []\n")
    monkeypatch.chdir(tmp_path)

    packaged = Path(server.__file__).parent / "presets.yaml"
    assert not packaged.exists()
    assert server._resolve_presets_path(packaged) == deployed


def test_explicit_presets_do_not_fall_back_to_deployment(monkeypatch, tmp_path):
    (tmp_path / "presets.yaml").write_text("presets: []\n")
    monkeypatch.chdir(tmp_path)

    explicit = tmp_path / "missing.yaml"
    assert server._resolve_presets_path(explicit) == explicit


def test_packaged_example_is_default_without_deployment(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)

    packaged = Path(server.__file__).parent / "presets.yaml"
    assert server._resolve_presets_path(packaged) == Path(str(packaged) + ".example")


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
