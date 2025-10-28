import logging
import pytest
import pandas as pd

from commec.control_list.containers import AccessionFormat, derive_accession_type, RegulationList, ListMode
import commec.control_list.list_data as ld
from commec.control_list.control_list import get_regulation, post_process_regulation_data

from commec.utils.logger import setup_console_logging

@pytest.mark.parametrize("new_list,expected_outcome", [
    pytest.param(*case) for case in [
        (RegulationList("list01","L1","www.l1.com",["NZ"],1), True), # Duplicate
        (RegulationList("list01","L2","www.l1.com",["NZ"],0), True), # Non-Duplicate
        (RegulationList("list01","L1","www.l2.com",["NZ"],0), False), # Wrong URL
        (RegulationList("list01","L1","www.l1.com",["NZ", "AU"],0), False), # Wrong Regions
        (RegulationList("list01","L1","www.l1.com",["AU"],0), False), # Wrong Region
        (RegulationList("list02","L1","www.l1.com",["NZ"],0), False), # Wrong Name
    ]
])
def test_list_overwrite_protections(new_list, expected_outcome):
    """
    Tests what occurs when two lists share an acronym that should not be concatenated.
    """
    setup_console_logging(logging.DEBUG)
    existing_list = RegulationList("list01","L1","www.l1.com",["NZ"],0)
    ld.add_regulated_list(existing_list)
    assert expected_outcome == ld.add_regulated_list(new_list)

def test_multiple_entrys():
    """
    Outcomes for where multiple lists contain duplicate taxid entries.
    """
    ld.clear()
    setup_console_logging(logging.DEBUG)
    assert ld.add_regulated_list(RegulationList("List01","L1","www.list1.com",["NZ"],ListMode.COMPLIANCE))
    assert ld.add_regulated_list(RegulationList("List02","L2","www.list2.com",["AU"],ListMode.COMPLIANCE))

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
    ld.add_regulated_taxid_data(input_list1_data)
    ld.add_regulated_taxid_data(input_list2_data)

    print("Preprocessing:")
    print(ld.REGULATED_TAXID_ANNOTATIONS.to_string())
    post_process_regulation_data()
    print("Postprocessing:")
    print(ld.REGULATED_TAXID_ANNOTATIONS.to_string())
    print(ld.REGULATED_TAXID_ANNOTATIONS.shape)

    # Twelve headings, however should only have 3 entries from the above 5 entries.
    assert ld.REGULATED_TAXID_ANNOTATIONS.shape == (3,12), "Incorrect number of imported Regulation Annotations."

    output = get_regulation(11320, AccessionFormat.TAXID)
    assert len(output) == 2, "Incorrect number of returned Regulations."

@pytest.mark.parametrize("accessions,expected_outcome", [
    pytest.param(*case) for case in [
        (["111","4","11084", 444], AccessionFormat.TAXID),
        (["HG992755.1","CAG2243592.1", "CP001814.1", "DI192294.1"], AccessionFormat.GENBANK),
        (["Q5VW38", "Q7L1I2", "V5XZS6", "B1P1E1"], AccessionFormat.UNIPROT),
    ]
])
def test_accession_identification_format(accessions, expected_outcome):
    """
    Tests some TaxIDs, Genbank records, and Uniprot accessions to ensure
    that they are correctly identified as such.
    Troublesome records should be added to the above list for continuous testing.
    """
    setup_console_logging(logging.DEBUG)
    for accession in accessions:
        outcome = derive_accession_type(accession)
        assert expected_outcome == outcome, f"{accession} failed to be identified. Expected {expected_outcome}, got {outcome}"