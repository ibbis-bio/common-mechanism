import logging
import pytest
import pandas as pd

import commec.regulation.containers as crc
from commec.regulation.regulation import get_regulation, post_process_regulation_data

from commec.utils.logger import setup_console_logging

@pytest.mark.parametrize("new_list,expected_outcome", [
    pytest.param(*case) for case in [
        (crc.RegulationList("list01","L1","www.l1.com",["NZ"],1), True), # Duplicate
        (crc.RegulationList("list01","L2","www.l1.com",["NZ"],0), True), # Non-Duplicate
        (crc.RegulationList("list01","L1","www.l2.com",["NZ"],0), False), # Wrong URL
        (crc.RegulationList("list01","L1","www.l1.com",["NZ", "AU"],0), False), # Wrong Regions
        (crc.RegulationList("list01","L1","www.l1.com",["AU"],0), False), # Wrong Region
        (crc.RegulationList("list02","L1","www.l1.com",["NZ"],0), False), # Wrong Name
    ]
])
def test_list_overwrite_protections(new_list, expected_outcome):
    """
    Tests what occurs when two lists share an acronym that should not be concatenated.
    """
    setup_console_logging(logging.DEBUG)
    existing_list = crc.RegulationList("list01","L1","www.l1.com",["NZ"],0)
    crc.add_regulated_list(existing_list)
    assert expected_outcome == crc.add_regulated_list(new_list)

def test_multiple_entrys():
    """
    Outcomes for where multiple lists contain duplicate taxid entries.
    """
    crc.clear()
    setup_console_logging(logging.DEBUG)
    assert crc.add_regulated_list(crc.RegulationList("List01","L1","www.list1.com",["NZ"],0))
    assert crc.add_regulated_list(crc.RegulationList("List02","L2","www.list2.com",["AU"],0))

    input_list1_data = pd.DataFrame([
        {
            "category": "Virus",
            "name": "Influenza A virus",
            "notes": "Seasonal pathogen",
            "derived_from": "",
            "preferred_taxonomy_name": "Influenza A virus",
            "other_taxonomy_name": "Flu A",
            "tax_id": "11320",
            "list_acronym": "L1",
            "target": "Human",
            "hazard_group": "HG2",
        },
        {
            "category": "Virus",
            "name": "Influenza B virus",
            "notes": "Seasonal pathogen",
            "derived_from": "",
            "preferred_taxonomy_name": "Influenza B virus",
            "other_taxonomy_name": "Flu B",
            "tax_id": "11321",
            "list_acronym": "L1",
            "target": "Human",
            "hazard_group": "HG2",
        },
    ])

    input_list2_data = pd.DataFrame([
        {
            "category": "Virus",
            "name": "Influenza A virus",
            "notes": "Seasonal pathogen",
            "derived_from": "",
            "preferred_taxonomy_name": "Influenza A virus",
            "other_taxonomy_name": "Flu A",
            "tax_id": "11320",
            "list_acronym": "L2",
            "target": "Human",
            "hazard_group": "HG2",
        },
        {   # Same as above, for L1, but differing input values, bare minimum
            #"name": "Influenza B virus",
            "tax_id": "11321",
            "list_acronym": "L1",
        },
        {   # Fully Duplicate entry
            "category": "Virus",
            "name": "Influenza B virus",
            "notes": "Seasonal pathogen",
            "derived_from": "",
            "preferred_taxonomy_name": "Influenza B virus",
            "other_taxonomy_name": "Flu B",
            "tax_id": "11321",
            "list_acronym": "L1",
            "target": "Human",
            "hazard_group": "HG2",
        },
    ])

    print("Input data 1")
    print(input_list1_data.to_string())
    print("Input data 2")
    print(input_list2_data.to_string())

    print("adding data...")
    crc.add_regulated_taxid_data(input_list1_data)
    crc.add_regulated_taxid_data(input_list2_data)

    print("Preprocessing:")
    print(crc.REGULATED_TAXID_ANNOTATIONS.to_string())
    post_process_regulation_data()
    print("Postprocessing:")
    print(crc.REGULATED_TAXID_ANNOTATIONS.to_string())
    print(crc.REGULATED_TAXID_ANNOTATIONS.shape)

    # Twelve headings, however should only have 3 entries from the above 5 entries.
    assert crc.REGULATED_TAXID_ANNOTATIONS.shape == (3,12), "Incorrect number of imported Regulation Annotations."

    output = get_regulation(11320, crc.AccessionFormat.TAXID)
    assert len(output) == 2, "Incorrect number of returned Regulations."