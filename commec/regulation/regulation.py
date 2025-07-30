# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

"""


"""
import os
from commec.regulation.containers import TaxidRegulation, RegulationList
from commec.regulation.initialisation import import_regulation_list



def clear(target : str | None = None) -> bool:
    """
    Removes the targeted list from the module state, or
    if no target is provided, clears the entired module state.
    returns whether operation was successful.
    """
    if target and target in _REG_LISTS:
        # implement targetted rmeoval logic.
        del _REG_LISTS[target]
        return True

    if not target:
        _REG_LISTS = None
        return True

    return False

def load(import_path : str | os.PathLike) -> bool:
    """
    Load a filepath recursively. Searches for valid "list folders"
    One found, the list folder is loaded into the modules state, and will
    """
    
    # check import_path for if info.txt exists.
    # Check if 

def is_regulated(taxid : int) -> bool:
