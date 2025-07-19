#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Static functions useful for dealing with common dictionary tasks.
"""

@staticmethod
def deep_update(to_update: dict[str, any], 
                has_updates: dict[str, any]) -> tuple[
                    dict[str,any], 
                    list[tuple[str,any]]]:
    """
    Recursively update a nested dictionary without completely overwriting nested dictionaries.
    Only already existing keys are updated. Any keys not existing in the dictionary
    to be updated are returned as a list of rejected key value pairs.
    -----
    Inputs:
    * to_update : dict[str, any] Dictionary to be updated.
    * has_updates : dict[str, any] New dictionary information to be added.
    ----
    Outputs:
    * updated : dict[str, any] a copy of the to_update dictionary, with values
    from any matching keys overridden by has_updates.
    * rejected : list[tuple[str,any] A list of the rejected key value pairs, i.e. 
    keys present in has_updates, but not present in to_update.
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
