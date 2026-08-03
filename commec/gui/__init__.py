# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""Web GUI subcommand for commec: a small Flask front-end for `commec screen`.

The subcommand contract (DESCRIPTION / add_args / run) lives in
`commec.gui.server`. This package initializer intentionally does not import
server, avoiding a double-import warning when it is run via
`python -m commec.gui.server`.
"""
