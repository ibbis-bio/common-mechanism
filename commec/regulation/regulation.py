# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

"""


"""
import os
import logging
from commec.regulation.containers import TaxidRegulation, RegulationList
import commec.regulation.containers as data
from commec.regulation.initialisation import import_regulations

logger = logging.getLogger(__name__)

def load(import_path : str | os.PathLike):
    """
    Load a filepath recursively. Searches for valid "list folders"
    One found, the list folder is loaded into the modules state, and will
    """

    if import_regulations(import_path):
        return
    
    # check for existance of sub folders and recurse on any present.
    for entry in os.scandir(import_path):
        if os.path.isdir(entry):
            load(entry)
    return

def get_regulation(taxid : int) -> list[tuple[RegulationList, TaxidRegulation]]:
    """
    Check the given TaxID against all imported regulated lists.
    The input taxid is 
    Any parent TaxIDS will also be recursively checked across 
    all regulated lists. The output is a list of every regulation 
    attributed to the original taxid and its parents, in the form of 
    a tuple containing the list, as well as the 
    taxid specific regulation information.
    """
    output_data : list[tuple[RegulationList, TaxidRegulation]] = []

    taxids_to_check = [taxid]
    taxid_parents_to_check = data.CHILD_TAXID_MAP[data.CHILD_TAXID_MAP["TaxID"] == taxid]["ParentTaxID"].to_list()
    taxids_to_check.append(taxid_parents_to_check)

    filtered_regulated_taxid_annotations = data.REGULATED_TAXID_ANNOTATIONS[
        data.REGULATED_TAXID_ANNOTATIONS["taxid"] in taxids_to_check]

    for _, row in filtered_regulated_taxid_annotations.iterrows():
        taxid_regulation_info = TaxidRegulation(
            row["taxonomy_category"],
            row["taxonomy_name"],
            row["notes"],
            row["preferred_taxonomy_name"],
            row["Taxid"],
            row["list_acronym"],
            row["target"],
            row["hazard_group"]
        )
        list_data = data.REGULATION_LISTS[taxid_regulation_info.list_acronym]
        output_data.append((list_data, taxid_regulation_info))

    return output_data


