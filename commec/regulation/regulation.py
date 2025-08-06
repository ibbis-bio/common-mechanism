# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

"""


"""
import os
import logging
from commec.regulation.containers import TaxidRegulation, RegulationList, REG_TAXID_LISTS, REG_LISTS
from commec.regulation.initialisation import import_regulations

logger = logging.getLogger(__name__)

def load(import_path : str | os.PathLike) -> bool:
    """
    Load a filepath recursively. Searches for valid "list folders"
    One found, the list folder is loaded into the modules state, and will
    """
    
    # check import_path for if info.txt exists.
    # if not, check for existance of sub folders and recurse.

    # if info.txt, import the file as a regulatedList object.
    # Check for existance of taxid info file, and import into List of RegulatedTaxids

    # inject both into the REG_TAXID_LISTS, and REG_LISTS data structures.

def is_regulated(taxid : int) -> bool:
    """
    Check the given TaxID against all imported regulated lists.
    Any parent TaxIDS will also be recursively checked across 
    all regulated lists. The output is a list of every regulation 
    attributed to the original taxid and its parents, in the form of 
    a tuple containing the list, as well as the 
    taxid specific regulation information.
    """
    output_data : list[tuple[RegulationList, TaxidRegulation]] = []

    for list_acronym, taxid_data in REG_TAXID_LISTS.items():
        taxid_regulation_data = taxid_data.get(taxid)
        if taxid_regulation_data:
            list_data = REG_LISTS.get(list_acronym)
            output_data.append((list_data, taxid_regulation_data))

            # Recursively check parents for their presence in all lists.
            if taxid_regulation_data.parent_taxid:
                parent_data = is_regulated(taxid_regulation_data.parent_taxid)
                if len(parent_data) > 0:
                    output_data.extend(parent_data)

    return output_data
