#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Utilities to set up commec package logging.
"""

import logging
import sys
import textwrap


class TextWrapFormatter(logging.Formatter):
    """
    Format multi-line log messages with proper vertical alignment, configurable styling,
    and text wrapping for longer messages.
    """

    def __init__(self, fmt=None, *args, continuation_marker="│ ", line_width=120, **kwargs):
        if fmt is None:
            fmt = f"%(levelname)-8s{continuation_marker}%(message)s"
        super().__init__(fmt, *args, **kwargs)
        self.continuation_marker = continuation_marker
        self.line_width = line_width

        # String to prepended to all lines of wrapped output except the first
        indent_size = self._find_message_start() - len(self.continuation_marker)
        self.indent = " " * indent_size + self.continuation_marker

    def _find_message_start(self):
        """
        Deterine how far to indent messages by formatting a dummy message.
        """
        sample = logging.LogRecord(
            name="dummy",
            level=logging.INFO,
            pathname="./test",
            lineno=0,
            msg="DUMMY_MESSAGE",
            args=(),
            exc_info=None,
        )
        sample.asctime = self.formatTime(sample)
        sample_formatted = super().format(sample)
        return sample_formatted.find(sample.msg)

    def format(self, record):
        message = super().format(record)
        lines = message.splitlines()

        formatted_lines = []
        # First line gets the levelname/timestamp/etc from super().format, then
        # long lines are wrapped with the indent
        wrapped_first = textwrap.wrap(
            lines[0],
            width=self.line_width,
            subsequent_indent=self.indent,
            break_long_words=False,
            break_on_hyphens=False,
        )
        formatted_lines.extend(wrapped_first)

        # When a message has newlines, lines after the first should be indented even if short
        for line in lines[1:]:
            wrapped = textwrap.wrap(
                line,
                width=self.line_width,
                initial_indent=self.indent,
                subsequent_indent=self.indent,
                break_long_words=False,
                break_on_hyphens=False,
            )
            formatted_lines.extend(wrapped)

        return "\n".join(formatted_lines)


def setup_console_logging(log_level=logging.INFO):
    """Set up logging to console."""
    commec_logger = logging.getLogger("commec")
    commec_logger.setLevel(log_level)

    # Check if the handler already exists to avoid duplicates
    if not any(isinstance(h, logging.StreamHandler) for h in commec_logger.handlers):
        console_handler = logging.StreamHandler()
        console_handler.setLevel(log_level)
        console_handler.setFormatter(TextWrapFormatter())
        commec_logger.addHandler(console_handler)

    add_logging_to_excepthook()


def setup_file_logging(filename, log_level=logging.INFO, log_mode="w"):
    """Set up logging to a file. Format determined based on level."""
    commec_logger = logging.getLogger("commec")

    # Ensure the logger level is set to the lowest level of any handler
    current_level = commec_logger.level or logging.INFO
    commec_logger.setLevel(min(current_level, log_level))

    # Log format has more detail if logging down to the debug level
    if log_level == logging.DEBUG:
        formatter = TextWrapFormatter(
            fmt="%(asctime)s│ %(levelname)-8s│ %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",  # Full ISO-like format
            line_width = 300, # Longer lines for debug purposes.
        )
    else:
        formatter = TextWrapFormatter("%(levelname)-8s│ %(message)s")

    # Update existing filehandlers, avoiding duplicates
    file_handler = None
    for handler in commec_logger.handlers:
        if (
            isinstance(handler, logging.FileHandler)
            and getattr(handler, "baseFilename", None) == filename
        ):
            file_handler = handler
            break

    file_handler = file_handler or logging.FileHandler(filename, log_mode)
    file_handler.setLevel(log_level)
    file_handler.setFormatter(formatter)
    commec_logger.addHandler(file_handler)


def add_logging_to_excepthook():
    """
    Ensure unhandled exceptions are logged to the commec package logger;
    original excepthook is still called.
    """
    original_excepthook = sys.excepthook

    def commec_exception_logger(exc_type, exc_value, exc_traceback):
        """Log exception to package logger."""
        commec_logger = logging.getLogger("commec")

        if commec_logger.handlers:
            # Log the exception message at ERROR level
            error_message = f"Unhandled exception: {exc_type.__name__}: {exc_value}"
            commec_logger.error(error_message)

            # Log the full traceback at the DEBUG level
            commec_logger.debug(
                "Exception traceback:", exc_info=(exc_type, exc_value, exc_traceback)
            )

        # Still call the original handler for console output
        original_excepthook(exc_type, exc_value, exc_traceback)

    sys.excepthook = commec_exception_logger


def set_log_level(log_level, update_only_handler_type=None):
    """
    Update the log level for the commec logger, as well as associated handlers.
    Optionally, restrict updates to only a particular class of handlers (e.g. StreamHandler).
    """
    commec_logger = logging.getLogger("commec")
    commec_logger.setLevel(log_level)

    handlers_to_update = commec_logger.handlers
    if update_only_handler_type:
        handlers_to_update = [
            h for h in handlers_to_update if isinstance(h, update_only_handler_type)
        ]

    for handler in handlers_to_update:
        handler.setLevel(log_level)
