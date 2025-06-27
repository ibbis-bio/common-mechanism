#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Static functions useful for dealing with common dictionary tasks.
"""
import logging
logger = logging.getLogger(__name__)

@staticmethod
def deep_update(to_update: dict, has_updates: dict):
    """
    Recursively update a nested dictionary without completely overwriting nested dictionaries.
    """
    rejected = []
    updated = to_update.copy()
    for key, value in has_updates.items():
        # If both values are dictionaries, recursively update
        if key in updated and isinstance(updated[key], dict) and isinstance(value, dict):
            updated[key], additional_rejects = deep_update(updated[key], value)
            rejected.extend(additional_rejects)
        # If not a dictionary, just copy the value.
        elif key in updated:
            updated[key] = value
        # If not present, we log an unexpected input one.
        else:
            rejected.append((key, value))
    return updated, rejected
