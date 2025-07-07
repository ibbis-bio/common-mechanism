#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
'''
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
'''

# Consider whether this can get away with being part of config. rename to IO config?

import json
import string
import os
from dataclasses import asdict, fields, is_dataclass
from typing import Dict, Type, get_origin, Any, get_args
from enum import StrEnum
from commec.config.result import ScreenResult, JSON_COMMEC_FORMAT_VERSION

class IoVersionError(RuntimeError):
    """Custom exception when handling differing versions with Commec output JSON."""

def encode_screen_data_to_json(input_screendata: ScreenResult,
                               output_json_filepath: string = "output.json") -> None:
    ''' Converts a ScreenResult class object into a JSON file at the given filepath.'''
    try:
        with open(output_json_filepath, "w", encoding="utf-8") as json_file:
            json.dump(asdict(input_screendata), json_file, indent=4)
    except TypeError as e:
        print("Error outputting JSON:", e)
        print(input_screendata)

def encode_dict_to_screen_data(input_dict : dict) -> ScreenResult:
    ''' Converts a dictionary into a ScreenResult object,
    any keys within the dictionary not part of the ScreenResult format are lost.
    any missing information will be simple set as defaults.'''
    return dict_to_dataclass(ScreenResult, input_dict)

# Convert the dictionary back to the dataclass or list of dataclass
def dict_to_dataclass(cls: Type, data: Dict[str, Any]) -> Any:
    '''
    Convert a dict, into appropriate dataclass, or list of dataclass,
    invalid keys to the dataclass structure are ignored.
    '''
    # Prepare a dictionary for filtered data
    filtered_data = {}

    if data is None:
        return filtered_data

    for f in fields(cls):
        field_name = f.name
        field_type = f.type

        if field_name in data:
            field_value = data[field_name]

            # Check if the field is a dataclass
            if is_dataclass(field_type):
                filtered_data[field_name] = dict_to_dataclass(field_type, field_value)
                continue

            # Check if the field is a list
            if get_origin(field_type) is list:
                item_type = get_args(field_type)[0]

                # Handle lists of StrEnums
                if issubclass(item_type, StrEnum):
                    filtered_data[field_name] = [item_type(item) for item in field_value]

                #Handles Dataclasses
                if is_dataclass(item_type) and isinstance(field_value, list):
                    filtered_data[field_name] = [
                        dict_to_dataclass(item_type, item) for item in field_value
                            if isinstance(item, dict)
                            and any(key in {f.name for f in fields(item_type)}
                            for key in item.keys()) or isinstance(item, item_type)]
                    continue

                filtered_data[field_name] = field_value
                continue

            # Check if the field is a dict of dataclasses
            if get_origin(field_type) is dict:
                _key_type, value_type = get_args(field_type)

                # Handle dicts of dataclasses
                if is_dataclass(value_type):
                    filtered_data[field_name] = {
                        key: dict_to_dataclass(value_type, value) if isinstance(value, dict)
                        else value for key, value in field_value.items()
                        if isinstance(value, (dict, value_type))
                    }
                    continue

                filtered_data[field_name] = field_value
                continue

            # Handle custom StrEnums
            if issubclass(field_type, StrEnum):
                try:
                    filtered_data[field_name] = field_type(field_value)
                except ValueError:
                    print(f"Invalid value '{field_value}' for "
                          f"field '{field_name}' of type {field_type}.")
                continue

            # Handle other field types
            filtered_data[field_name] = field_value

    # Create an instance of the dataclass with the filtered data
    return cls(**filtered_data)

def get_screen_data_from_json(input_json_filepath: string) -> ScreenResult:
    ''' Loads a JSON file from given filepath and returns
    a populated ScreenResult object from its contents. If the file does not
    exist, then returns a new screen data object.'''
    if not os.path.exists(input_json_filepath):
        return ScreenResult()

    json_string : str
    with open(input_json_filepath, "r", encoding="utf-8") as json_file:
        # Read the file contents as a string
        json_string = json_file.read()
    my_data : dict = json.loads(json_string)

    # Check version of imported json.
    input_version = my_data["commec_info"]["json_output_version"]
    if not input_version == JSON_COMMEC_FORMAT_VERSION:
        raise IoVersionError(f"Version difference between input (v.{input_version}) and"
                            f" expected (v.{JSON_COMMEC_FORMAT_VERSION})"
                            f": {input_json_filepath}")
    return encode_dict_to_screen_data(my_data)
