import pytest
import os
from commec.json_io import *
from dataclasses import asdict

@pytest.fixture
def test_screendata():
    '''Fixture to provide the ScreenData for testing.'''
    return ScreenData(
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

@pytest.fixture
def empty_screendata():
    '''Fixture to provide the ScreenData for testing.'''
    return ScreenData(recommendation="PASS")

@pytest.mark.parametrize("test_data_fixture",["test_screendata", "empty_screendata"])
def test_json_io(tmp_path, request, test_data_fixture):
    ''' Test to ensure that read/write for JSON ScreenData I/O is working correctly.'''
    test_data = request.getfixturevalue(test_data_fixture)
    json_filename1 = tmp_path / "testread1.json"
    json_filename2 = tmp_path / "testread2.json"
    encode_screen_data_to_json(test_data, json_filename1)
    test_data_retrieved = get_screen_data_from_json(json_filename1)
    encode_screen_data_to_json(test_data_retrieved, json_filename2)
    test_data_retrieved_twice = get_screen_data_from_json(json_filename2)

    # Convert both original and retrieved data to dictionaries and compare
    assert asdict(test_data) == asdict(test_data_retrieved), (
        f"JSON Write/Read interpreter failed.\n"
        f"Test JSON Reference data: \n{asdict(test_data)}\n"
        f"Test JSON output data: \n{asdict(test_data_retrieved)}"
    )

    # Convert both original and retrieved data to dictionaries and compare
    assert asdict(test_data) == asdict(test_data_retrieved_twice), (
        f"JSON Write/Read/Write/Read interpreter failed.\n"
        f"Test JSON Reference data: \n{asdict(test_data)}\n"
        f"Test JSON output data: \n{asdict(test_data_retrieved)}"
    )


def test_erroneous_info(tmp_path, test_screendata):
    ''' Test to ensure that read/write for JSON ScreenData I/O is working correctly.'''
    test_data = test_screendata
    json_filename3 = tmp_path / "testread3.json"
    json_filename4 = tmp_path / "testread4.json"

    encode_screen_data_to_json(test_data, json_filename3)
    test_data_retrieved = get_screen_data_from_json(json_filename3)

    # Add erroneous information
    test_data_dict = asdict(test_data_retrieved)
    test_data_dict["ExtraStuff1"] = "ExtraBitStuff1"
    test_data_dict["biorisks"][0]["ExtraStuff2"] = "ExtraBitStuff2"
    test_data_dict["biorisks"].append("ExtraStuff3")
    test_data_dict["biorisks"].append({"ExtraDictStuff4" : 9999})
    test_data_dict2 = encode_dict_to_screen_data(test_data_dict)
    encode_screen_data_to_json(test_data_dict2, json_filename4)
    test_data_retrieved = get_screen_data_from_json(json_filename4)

    # Convert both original and retrieved data to dictionaries and compare
    assert asdict(test_data) == asdict(test_data_retrieved), (
        f"JSON Write/Read interpreter failed.\n"
        f"Test JSON Reference data: \n{asdict(test_data)}\n"
        f"Test JSON output data: \n{asdict(test_data_retrieved)}"
    )
