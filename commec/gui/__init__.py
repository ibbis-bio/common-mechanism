# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""Web GUI subcommand for commec: a small Flask front-end for `commec screen`.

The subcommand contract (DESCRIPTION / add_args / run) lives in
`commec.gui.server`. commec.cli imports it defensively (guarded by try/except)
so a missing Flask disables only `commec gui` rather than the whole CLI. This
package initializer intentionally does NOT import server: it keeps
`import commec.gui` Flask-free and avoids a double-import warning when server is
run via `python -m commec.gui.server`.
"""
