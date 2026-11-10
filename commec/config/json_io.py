#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Set of tools for retrieving and storing information important to screen
 outputs. Information is stored as a structure of dataclasses (headed by ScreenResult),
 and are converted between the dataclass / dict / json_file as required. The
 conversions are done dynamically, and it is recommended to only use and
 interact with the dataclasses only, to maintain version format, and not
 create erroneous outputs to the JSON which wont be read back in. This
 ensures an expected i/o behaviour.

 The single exception to this is the "annotations" dictionary, present
 in the HitDescription, which contains non-structured information, and is
 populated with differing information under differing keys depending on
 which step the information is derived (Biorisk, Taxonomy etc)

 In this way, the JSON object serves as a common state, that can be updated
 whilst not being temporally appended like a log file i.e. .screen file.

 The JSON stores all pertinent information of a run.
"""

# Consider whether this can get away with being part of config. rename to IO config?

import json
import os
import string
from dataclasses import asdict, fields, is_dataclass
from enum import StrEnum
from types import UnionType
from typing import Any, Dict, Type, get_args, get_origin

from commec.config.result import (
    ASSEMBLY_JSON_FORMAT_VERSION,
    JSON_COMMEC_FORMAT_VERSION,
    ScreenResult,
)


class IoVersionError(RuntimeError):
    """Custom exception when handling differing versions with Commec output JSON."""


def encode_screen_data_to_json(
    input_result: ScreenResult, output_json_filepath: string = "output.json"
) -> None:
    """Converts a ScreenResult class object into a JSON file at the given filepath."""
    try:
        document = asdict(input_result)
        if input_result.assembly_info is None:
            document.pop("assembly_info", None)
            document.pop("assembly_summary", None)
            for query in document["queries"].values():
                query.pop("assembly", None)
                query.pop("combined_status", None)
        with open(output_json_filepath, "w", encoding="utf-8") as json_file:
            json.dump(document, json_file, indent=2)
    except TypeError as e:
        print("Error outputting JSON:", e)
        print(input_result)


def encode_dict_to_screen_data(input_dict: dict) -> ScreenResult:
    """Converts a dictionary into a ScreenResult object,
    any keys within the dictionary not part of the ScreenResult format are lost.
    any missing information will be simple set as defaults."""
    return dict_to_dataclass(ScreenResult, input_dict)


# Convert the dictionary back to the dataclass or list of dataclass
def dict_to_dataclass(cls: Type, data: Dict[str, Any]) -> Any:
    """
    Convert a dict, into appropriate dataclass, or list of dataclass,
    invalid keys to the dataclass structure are ignored.
    """
    # Prepare a dictionary for filtered data
    filtered_data = {}

    if data is None:
        return None
    if get_origin(cls) in (UnionType,):
        non_null = [item for item in get_args(cls) if item is not type(None)]
        cls = non_null[0] if len(non_null) == 1 else cls

    for f in fields(cls):
        field_name = f.name
        field_type = f.type

        if field_name in data:
            field_value = data[field_name]

            # Unwrap optional dataclasses before checking the dataclass itself.
            if get_origin(field_type) is UnionType:
                non_null = [
                    item for item in get_args(field_type) if item is not type(None)
                ]
                if len(non_null) == 1 and is_dataclass(non_null[0]):
                    filtered_data[field_name] = (
                        dict_to_dataclass(non_null[0], field_value)
                        if field_value is not None
                        else None
                    )
                    continue

            # Check if the field is a dataclass
            if is_dataclass(field_type):
                filtered_data[field_name] = dict_to_dataclass(field_type, field_value)
                continue

            # Check if the field is a list
            if get_origin(field_type) is list:
                item_type = get_args(field_type)[0]

                # Handle lists of StrEnums
                if issubclass(item_type, StrEnum):
                    filtered_data[field_name] = [
                        item_type(item) for item in field_value
                    ]

                # Handles Dataclasses
                if is_dataclass(item_type) and isinstance(field_value, list):
                    filtered_data[field_name] = [
                        dict_to_dataclass(item_type, item)
                        for item in field_value
                        if isinstance(item, dict)
                        and any(
                            key in {f.name for f in fields(item_type)}
                            for key in item.keys()
                        )
                        or isinstance(item, item_type)
                    ]
                    continue

                filtered_data[field_name] = field_value
                continue

            # Check if the field is a dict of dataclasses
            if get_origin(field_type) is dict:
                _key_type, value_type = get_args(field_type)

                # Handle dicts of dataclasses
                if is_dataclass(value_type):
                    filtered_data[field_name] = {
                        key: dict_to_dataclass(value_type, value)
                        if isinstance(value, dict)
                        else value
                        for key, value in field_value.items()
                        if isinstance(value, (dict, value_type))
                    }
                    continue

                filtered_data[field_name] = field_value
                continue

            if get_origin(field_type) is UnionType:
                if field_value is None:
                    filtered_data[field_name] = None
                    continue
                non_null = [
                    item for item in get_args(field_type) if item is not type(None)
                ]
                field_type = non_null[0] if len(non_null) == 1 else field_type

            # Handle custom StrEnums
            if isinstance(field_type, type) and issubclass(field_type, StrEnum):
                try:
                    filtered_data[field_name] = field_type(field_value)
                except ValueError:
                    print(
                        f"Invalid value '{field_value}' for "
                        f"field '{field_name}' of type {field_type}."
                    )
                continue

            # Handle other field types, including optional enums.
            if get_origin(field_type) is UnionType:
                non_null = [
                    item for item in get_args(field_type) if item is not type(None)
                ]
                if (
                    len(non_null) == 1
                    and isinstance(non_null[0], type)
                    and issubclass(non_null[0], StrEnum)
                ):
                    filtered_data[field_name] = non_null[0](field_value)
                    continue
            filtered_data[field_name] = field_value

    # Create an instance of the dataclass with the filtered data
    return cls(**filtered_data)


def get_screen_data_from_json(input_json_filepath: string) -> ScreenResult:
    """Loads a JSON file from given filepath and returns
    a populated ScreenResult object from its contents. If the file does not
    exist, then returns a new screen data object."""
    if not os.path.exists(input_json_filepath):
        return ScreenResult()

    json_string: str
    with open(input_json_filepath, "r", encoding="utf-8") as json_file:
        # Read the file contents as a string
        json_string = json_file.read()
    my_data: dict = json.loads(json_string)

    # Version 0.6 lacks assembly evidence. Dataclass defaults preserve its behaviour.
    input_version = my_data["commec_info"]["json_output_version"]
    if input_version not in {JSON_COMMEC_FORMAT_VERSION, ASSEMBLY_JSON_FORMAT_VERSION}:
        raise IoVersionError(
            f"Version difference between input (v.{input_version}) and"
            f" expected (v.{JSON_COMMEC_FORMAT_VERSION})"
            f": {input_json_filepath}"
        )
    return encode_dict_to_screen_data(my_data)
