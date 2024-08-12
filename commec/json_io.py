#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
'''
   Set of tools for retrieving and storing information important to screen
    outputs. Information is stored as a structure of dataclasses, and are 
    converted between the dataclass / dict / json_file as required. The 
    conversions are done dynamically, and it is recommended to only use and 
    interact with the dataclasses only, to maintain version format, and not 
    create erroneous outputs to the JSON which wont be read back in. This
    ensures an expected i/o behaviour.

   In this way, the JSON object serves as a common state, that can be updated
    whilst not being temporally appended like a log file i.e. .screen file.
'''

import json
import string
from dataclasses import dataclass, asdict, fields, field, is_dataclass
import os
from typing import Dict, Type, get_origin, Any, get_args

# Seperate versioning for the output JSON.
JSON_COMMEC_FORMAT_VERSION = "1.0"

@dataclass
class MatchRange:
    '''Container for information on how a match maps to a query, match direction
      is +1, or positive, if the query and match are in the same direction, and 
      -1 if they are inverse, and 0 if not set.'''
    match_start : int = 0
    match_end : int = 0
    query_start : int = 0
    query_end : int = 0
    query_direction : int = 0

@dataclass
class BenignData:
    '''Container to hold data related to the benign status of a taxonomy search'''
    benign : bool = False
    benign_match_type : str = ""
    percent_identity_e_value : float = 0.0

@dataclass
class MatchFields:
    '''Container to hold information for a match to the query identified in a commec screen. 
    molecule_alphabet is expected to contain either 'aa' for aminoacids, or 'nt' for nucelotides.
    range contains information pertaining to over what range the match occured.'''
    molecule_alphabet : str = ""
    description : str = ""
    id : int = 0
    taxon : str = ""
    kingdom : str = ""
    range : MatchRange = field(default_factory = MatchRange)
    is_benign : BenignData = field(default_factory = BenignData)
    # etc etc

@dataclass
class BioRiskData:
    '''Container dataclass for a list of matches to biorisks 
    identified from a commec database screen.'''
    matches : list[MatchFields] = field(default_factory = list)

@dataclass
class TaxonomyData:
    '''Container dataclass for a list of matches of taxonomy hits, 
    identified from a commec database screen.'''
    is_regulated : bool
    matches : list[MatchFields] = field(default_factory = list)

@dataclass
class QueryData:
    '''Container to hold data related to the Query
      used in a commec database screen'''
    name : str = ""
    length : int = 0
    sequence : str = ""

@dataclass
class CommecRunInformation:
    '''Container dataclass to hold general run information for a commec screen '''
    commec_version : str = "0.1.2"
    json_output_version : str = JSON_COMMEC_FORMAT_VERSION
    # add other settings / run parameters.

@dataclass
class ScreenData:
    ''' Root dataclass to hold all data related to the screening of an individual query by commec.'''
    recommendation : str = ""
    query : QueryData = field(default_factory = QueryData)
    biorisks : list[BioRiskData] = field(default_factory = list)
    taxonomies : list[TaxonomyData] = field(default_factory = list)
    commec_info : CommecRunInformation = field(default_factory = CommecRunInformation)

    def format(self):
        ''' Format this ScreenData as a json string to pass to a standard out if desired.'''
        return str(asdict(self))

# The above could be moved to a custom .py script for variable importing under version control.

def encode_screen_data_to_json(input_screendata: ScreenData, output_json_filepath: string = "output.json"):
    ''' Converts a ScreenData class object into a JSON file at the given filepath.'''
    with open(output_json_filepath, "w", encoding="utf-8") as json_file:
        json.dump(asdict(input_screendata), json_file, indent=4)

def encode_dict_to_screen_data(input_dict : dict):
    ''' Converts a dictionary into a ScreenData object, 
    any keys within the dictionary not part of the ScreenData format are lost.
    any missing information will be simple set as defaults.'''
    return dict_to_dataclass(ScreenData, input_dict)



# Convert the dictionary back to the dataclass or list of dataclass
def dict_to_dataclass(cls: Type, data: Dict[str, Any]) -> Any:
    ''' 
    Convert a dict, into appropriate dataclass, or list of dataclass, 
    invalid keys to the dataclass structure are ignored.
    '''
    # Prepare a dictionary for filtered data
    filtered_data = {}

    for f in fields(cls):
        field_name = f.name
        field_type = f.type
        if field_name in data:
            field_value = data[field_name]

            # Check if the field is a dataclass
            if is_dataclass(field_type):
                filtered_data[field_name] = dict_to_dataclass(field_type, field_value)
                continue

            # Check if the field is a list of dataclasses
            if get_origin(field_type) is list:
                item_type = get_args(field_type)[0]
                if is_dataclass(item_type) and isinstance(field_value, list):
                    filtered_data[field_name] = [
                        dict_to_dataclass(item_type, item) for item in field_value
                        if isinstance(item, dict) 
                            and any(key in {f.name for f in fields(item_type)} for key in item.keys())
                            or isinstance(item, item_type)]
                    continue
                filtered_data[field_name] = field_value
                continue

            # Handle other field types
            filtered_data[field_name] = field_value

    # Create an instance of the dataclass with the filtered data
    return cls(**filtered_data)

def get_screen_data_from_json(input_json_filepath: string):
    ''' Loads a JSON file from given filepath and returns a populated ScreenData object from its contents.'''
    json_string : str
    with open(input_json_filepath, "r", encoding="utf-8") as json_file:
        # Read the file contents as a string
        json_string = json_file.read()
    my_data : dict = json.loads(json_string)

    # Check version of imported json.
    input_version = my_data["commec_info"]["json_output_version"]
    if not input_version == JSON_COMMEC_FORMAT_VERSION:
        raise RuntimeError("""Version difference between input (v.{input_version}) and
                            expected (v.{JSON_COMMEC_FORMAT_VERSION}) JSON state file: {input_json_filepath}""")
    return encode_dict_to_screen_data(my_data)

def test1():
    ''' Simple test to ensure that read/write for json ScreenData I/O is working correctly.'''
    print("Testing JSON i/o.")

    test_data : ScreenData = ScreenData(
        recommendation="PASS",
        query = QueryData(
            name="Query1",
            length=10,
            sequence="ABCDEFGHIJ"
        ),
        biorisks = [
            BioRiskData(
                matches = [
                    MatchFields(
                        "nt"
                    )
                ]
            ),
            BioRiskData(

            ),
        ],
        taxonomies = [

        ],
        commec_info = CommecRunInformation(
            commec_version="0.1.2",
            json_output_version=JSON_COMMEC_FORMAT_VERSION,
        )
    )

    encode_screen_data_to_json(test_data, "test.json")
    test_data_retrieved = get_screen_data_from_json("test.json")

    test_dict_1 = str(asdict(test_data))
    test_dict_2 = str(asdict(test_data_retrieved))
    if test_dict_1 == test_dict_2:
        print("JSON interpreter T1 working correctly.")
        return 0
    print("JSON interpreter FAILED.")
    print("Test JSON Reference data: \n", test_dict_1)
    print("Test JSON output data: \n", test_dict_2)
    return 1

def test2():
    ''' Simple test to ensure that read/write for json ScreenData I/O is working correctly.'''
    test_data : ScreenData = ScreenData()
    test_data.recommendation = "PASS"
    encode_screen_data_to_json(test_data, "test.json")
    test_data_retrieved = get_screen_data_from_json("test.json")
    test_dict_1 : str = test_data.format()
    test_dict_2 : str = test_data_retrieved.format()
    if test_dict_1 == test_dict_2:
        print("JSON interpreter T2 working correctly.")
        return 0
    print("JSON interpreter FAILED.")
    print("Test JSON Reference data: \n", test_dict_1)
    print("Test JSON output data: \n", test_dict_2)
    return 1

def test3():
    ''' Simple test to ensure that read/write for json ScreenData I/O is working correctly.'''
    test_data : ScreenData = ScreenData()
    test_data.recommendation = "PASS"
    encode_screen_data_to_json(test_data, "test.json")
    test_data_retrieved = get_screen_data_from_json("test.json")
    test_dict_1 = str(test_data.recommendation)
    test_dict_2 = str(test_data_retrieved.recommendation)
    if test_dict_1 == test_dict_2:
        print("JSON interpreter T3 working correctly.")
        return 0

    print("JSON interpreter FAILED.")
    print("Test JSON Reference data: \n", test_dict_1)
    print("Test JSON output data: \n", test_dict_2)
    return 1

def test4():
    ''' Simple test to ensure that read/write for json ScreenData I/O is working correctly.'''
    test_data : ScreenData = ScreenData(
        recommendation="PASS",
        query = QueryData(
            name="Query1",
            length=10,
            sequence="ABCDEFGHIJ"
        ),
        biorisks = [
            BioRiskData(
                matches = [
                    MatchFields(
                        "nt"
                    )
                ]
            ),
        ],
        taxonomies = [

        ],
        commec_info = CommecRunInformation(
            commec_version="0.1.2",
            json_output_version=JSON_COMMEC_FORMAT_VERSION,
        )
    )

    encode_screen_data_to_json(test_data, "test.json")
    test_data_retrieved = get_screen_data_from_json("test.json")
    test_data_dict = asdict(test_data_retrieved)
    test_data_dict["ExtraStuff1"] = "ExtraBitStuff1"
    test_data_dict["biorisks"][0]["ExtraStuff2"] = "ExtraBitStuff2"
    test_data_dict["biorisks"].append("ExtraStuff3")
    test_data_dict["biorisks"].append({"ExtraDictStuff4" : 9999})
    test_data_dict2 = encode_dict_to_screen_data(test_data_dict)
    encode_screen_data_to_json(test_data_dict2, "test2.json")
    test_data_retrieved = get_screen_data_from_json("test2.json")

    test_dict_1 = str(asdict(test_data))
    test_dict_2 = str(asdict(test_data_dict2))
    if test_dict_1 == test_dict_2:
        print("JSON interpreter T4 working correctly.")
        return 0
    print("JSON interpreter FAILED.")
    print("Test JSON Reference data: \n", test_dict_1)
    print("Test JSON output data: \n", test_dict_2)
    return 1

if __name__ == "__main__":
    # Use main as unit test entry point.
    failed_tests : int = 0
    failed_tests += test1() # O/I test.
    failed_tests += test2() # Missing values test
    failed_tests += test3() # Missing values test
    failed_tests += test4() # Extra value test
    # Added values test
    # Version mismatch test
    print("Failed", failed_tests, "out of 4 tests.")
    if failed_tests == 0:
        print("Testing SUCCESSFUL")
        os.remove("test.json")
        os.remove("test2.json")
    else:
        print("Testing FAILED")

