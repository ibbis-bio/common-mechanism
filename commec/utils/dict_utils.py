#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Static functions useful for dealing with common dictionary tasks.
"""

@staticmethod
def deep_update(to_update: dict, has_updates: dict):
    """
    Recursively update a nested dictionary without completely overwriting nested dictionaries.
    """
    updated = to_update.copy()
    for key, value in has_updates.items():
        # If both values are dictionaries, recursively update
        if key in updated and isinstance(updated[key], dict) and isinstance(value, dict):
            updated[key] = deep_update(updated[key], value)
        else:
            updated[key] = value
    return updated
