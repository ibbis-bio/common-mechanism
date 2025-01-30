import pytest
from dataclasses import asdict
from commec.config.json_io import *
from commec.config.result import *
from commec.tools.search_handler import SearchToolVersion

@pytest.fixture
def test_screendata():
    '''Fixture to provide the ScreenResult for testing.'''
    return ScreenResult(
        #recommendation="PASS",
        commec_info = ScreenRunInfo(
            commec_version="0.1.2",
            json_output_version=JSON_COMMEC_FORMAT_VERSION,
            biorisk_database_info=SearchToolVersion("HMM 0.0.0","DB 0.0.0"),
            protein_database_info=SearchToolVersion("Blast 0.0.0","DB 0.0.0"),
            nucleotide_database_info=SearchToolVersion("Blast 0.0.0","DB 0.0.0"),
            benign_protein_database_info=SearchToolVersion("Blast 0.0.0","DB 0.0.0"),
            benign_rna_database_info=SearchToolVersion("Blast 0.0.0","DB 0.0.0"),
            benign_synbio_database_info=SearchToolVersion("Blast 0.0.0","DB 0.0.0"),
            time_taken="00:00:00:00",
            date_run="1.1.2024",
        ),
        queries= {
            "Query1":
            QueryResult(
                query_name="Query1",
                query_length=10,
                sequence="ABCDEFGHIJ",
                recommendation = QueryRecommendationContainer(),
                hits = {
                    "ImportantProtein1":
                    HitResult(
                        recommendation=HitRecommendationContainer(Recommendation.WARN, ScreenStep.BIORISK),
                        name="ImportantProtein1",
                        annotations = {"domain" : ["Bacteria"]},
                        ranges = [
                            MatchRange(
                                e_value = 0.0,
                                match_start = 0,
                                match_end = 10,
                                query_start = 0,
                                query_end = 10
                            )
                        ]
                    )
                }
            )
        },
    )

@pytest.fixture
def empty_screendata():
    '''Fixture to provide the ScreenResult for testing.'''
    return ScreenResult()

@pytest.mark.parametrize("test_data_fixture",["test_screendata", "empty_screendata"])
def test_json_io(tmp_path, request, test_data_fixture):
    ''' Test to ensure that read/write for JSON ScreenResult I/O is working correctly.'''
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
    ''' Test to ensure that read/write for JSON ScreenResult I/O is working correctly.'''
    test_data = test_screendata
    json_filename3 = tmp_path / "testread3.json"
    json_filename4 = tmp_path / "testread4.json"

    encode_screen_data_to_json(test_data, json_filename3)
    test_data_retrieved = get_screen_data_from_json(json_filename3)

    # Add erroneous information
    test_data_dict = asdict(test_data_retrieved)
    test_data_dict["ExtraStuff1"] = "ExtraBitStuff1"
    test_data_dict["queries"]["Query1"]["ExtraStuff2"] = "ExtraBitStuff2"
    test_data_dict["queries"]["Query1"]["hits"]["ImportantProtein1"]["ranges"].append("ExtraStuff3")
    test_data_dict["queries"]["Query1"]["hits"]["ImportantProtein1"]["ranges"].append({"ExtraDictStuff4" : 9999})
    test_data_dict2 = encode_dict_to_screen_data(test_data_dict)
    encode_screen_data_to_json(test_data_dict2, json_filename4)
    test_data_retrieved = get_screen_data_from_json(json_filename4)

    # Convert both original and retrieved data to dictionaries and compare
    assert asdict(test_data) == asdict(test_data_retrieved), (
        f"JSON Write/Read interpreter failed.\n"
        f"Test JSON Reference data: \n{asdict(test_data)}\n\n\n\n"
        f"Test JSON output data: \n{asdict(test_data_retrieved)}\n\n\n\n"
    )

def test_recommendation_ordering():
    assert Recommendation.PASS.importance < Recommendation.FLAG.importance
    assert compare(Recommendation.PASS, Recommendation.FLAG) == Recommendation.FLAG

def test_adding_data_to_existing():
    """
    Tests to ensure the mutability of writing to queries is working as expected.
    """
    def write_info(input_query : QueryResult):
        input_query.recommendation.biorisk_screen = Recommendation.PASS
    
    new_screen_data = ScreenResult()
    new_screen_data.queries["test01"] = QueryResult("test01", 10, "ATGCATGCAT", Recommendation.FLAG)
    write_query = new_screen_data.get_query("test01")
    write_info(write_query)
    assert new_screen_data.queries["test01"].recommendation.biorisk_screen == Recommendation.PASS
