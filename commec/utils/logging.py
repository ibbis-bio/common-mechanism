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
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(logging.Formatter("%(levelname)-8s | %(message)s"))
    commec_logger.addHandler(console_handler)

    add_logging_to_excepthook()


def setup_file_logging(self, filename, log_level=logging.INFO, log_mode="w"):
    """Set up logging to a file. Format determined based on level."""
    commec_logger = logging.getLogger("commec")
    file_handler = logging.FileHandler(filename, log_mode)
    file_handler.setLevel(log_level)

    formatter = logging.Formatter("%(levelname)-8s | %(message)s")
    # Add more info if logging down to the debug level
    if log_level == logging.DEBUG:
        formatter = logging.Formatter(
            "%(asctime)s | %(levelname)-8s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",  # Full ISO-like format
        )

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
