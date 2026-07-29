"""Input handling for the Commec GUI."""

from commec.gui import server


def test_short_sequences_are_passed_to_commec(tmp_path, monkeypatch):
    """The GUI leaves sequence-length handling to the screening pipeline."""
    monkeypatch.setitem(server.CFG, "password_hash", None)
    monkeypatch.setitem(server.CFG, "require_local_auth", False)
    monkeypatch.setitem(server.CFG, "runs_dir", tmp_path)
    monkeypatch.setitem(
        server.CFG,
        "presets",
        [{"id": "test", "label": "Test", "config": {}}],
    )
    monkeypatch.setattr(server, "_missing_databases", lambda _preset: [])
    monkeypatch.setattr(server.threading.Thread, "start", lambda _thread: None)
    server.JOBS.clear()
    server.app.config.update(TESTING=True)

    response = server.app.test_client().post(
        "/screen",
        data={
            "preset": "test",
            "label": "short-input",
            "sequence_text": ">short\nACGT\n",
            # Older clients may still submit the removed field. It must not
            # restore the GUI's former length filter.
            "skip_short": "1",
        },
        environ_overrides={"REMOTE_ADDR": "127.0.0.1"},
    )

    assert response.status_code == 200
    run_dir = next(path for path in tmp_path.iterdir() if path.is_dir())
    assert (run_dir / "input.fasta").read_text(encoding="utf-8") == ">short\nACGT\n"
    server.JOBS.clear()
