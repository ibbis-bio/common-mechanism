#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Utilities to set up commec package logging.
"""

import logging
import sys


def setup_console_logging(log_level=logging.INFO):
    """Set up logging to console."""
    commec_logger = logging.getLogger("commec")
    commec_logger.setLevel(log_level)

    # Check if the handler already exists to avoid duplicates
    if not any(isinstance(h, logging.StreamHandler) for h in commec_logger.handlers):
        console_handler = logging.StreamHandler()
        console_handler.setLevel(log_level)
        console_handler.setFormatter(logging.Formatter("%(levelname)-8s | %(message)s"))
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
        formatter = logging.Formatter(
            "%(asctime)s | %(levelname)-8s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",  # Full ISO-like format
        )
    else:
        formatter = logging.Formatter("%(levelname)-8s | %(message)s")

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
