# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

"""
Scripts that control the initialisation of the regulations
module state and imports regulation lists.

Import is handled recursively.

A valid list folder has the following layout:

regulated-lists/  		# This is the level at which pass to yaml / cli
├── uk-coshh/			# Arbitrary filename, contains 1x list.  
├── austgroup/          # The valid list folder layout:
│   ├── regulated_taxids.csv	# Taxid List annotations.
│   ├── regulated_taxids_and_children.csv
|   |		“
|   |		TaxID, ParentTaxID
|   |		”
│   ├── info.csv		# List of regions this list affects.
│   │   └── “
|   |		full_list_name,list_acronym,list_url,region_name, region_acronym
|   |		Special list, SL, www.sl.com, European Union, EU
|   |		”
│   └── ...
└── ...

"""
import os
import logging
import pandas as pd

from commec.regulation.containers import (
    RegulationList,
    Region
    )

import commec.regulation.containers as rc

logger = logging.getLogger(__name__)

def is_valid_regulation_list_folder(input_path : str | os.PathLike) -> bool:
    """
    Checks if the supplies folder is a valid regulation list.
    i.e. contains the regulated_taxids.csv, regions_info.csv
    """
    info_filename = os.path.join(input_path, "info.csv")
    data_filename = os.path.join(input_path, "regulated_taxids.csv")
    child_lut_filename = os.path.join(input_path, "children_of_regulated_taxids.csv")

    if not os.path.isdir(input_path): return False
    if not os.path.isfile(info_filename): return False
    if not os.path.isfile(data_filename): return False
    if not os.path.isfile(child_lut_filename): return False
    return True

def _import_regulation_list_info(input_path : str | os.PathLike):
    """
    Imports the info.csv file from a regulation list folder.
    Ensures that each regulation list information is appropriate keyed by
    acronym into the REG_LISTS data container.
    Ensures that existing regulation lists are not overwritten, and 
    warns the user if overwritting unique data (i.e. acroynym clash) has occured.
    """

    list_info = pd.read_csv(input_path)
    for _, row in list_info.iterrows():

        new_list = RegulationList(
            row["full_list_name"],
            row["list_acronym"],
            row["list_url"],
            list(Region(name = row["region_name"],
                    acronym = row["region_code"])))

        # Check list doesn't already exist, or is not overwritting another.
        list_key = row["list_acronym"]
        if rc.REGULATION_LISTS.get(list_key):
            existing_list = rc.REGULATION_LISTS[list_key]
            if new_list == existing_list:
                logger.debug("List already exists, no need to add.")
                continue

            # Check that these are different lists by comparing the names.
            logger.error("Error when importing list, shared acronym"
                         " between regulated list A and B, List B will be skipped.")
            continue

        # Add list.
        rc.REGULATION_LISTS[list_key] = new_list

def _import_regulation_taxid_data(input_path : str | os.PathLike):
    """
    Imports annotated regulated taxid information from the regulated_taxids.csv file
    within a regulated list provided to commec.
    Each regulated taxid is checked to ensure the regulated list is valid before
    inclusion.
    Concatenates the regulated taxid info into the global dataframe.
    """
    taxid_info = pd.read_csv(input_path)

    # Only include data whose list acronym exists.
    mask = taxid_info["list_acronym"].apply(
        lambda list_key: rc.REGULATION_LISTS.get(list_key) is not None)
    valid_list_taxid_info = taxid_info[mask]

    # Warn about dropped rows
    dropped = taxid_info[~mask]
    if not dropped.empty:
        logger.warning("The following list acronyms were not valid from %s", input_path)
        logger.warning(dropped["list_acronym"].tolist())

    # Append the new list data:
    rc.add_regulated_taxid_data(valid_list_taxid_info)

def _import_child_to_regulated_taxid_relationship(input_path : str | os.PathLike):
    """
    Imports child to regulated taxid look up data data from the 
    children_of_regulated_taxids.csv file within a regulated list provided to commec.
    Concatenates the child LUT info into the global dataframe.
    """
    child_lut = pd.read_csv(input_path)
    # Append the new list data:
    rc.add_regulated_taxid_data(child_lut)

def import_regulations(input_path : str | os.PathLike) -> bool:
    """
    Given a directory, checks that the directory is valid,
    If the directory is valid, it will load the regulation data.
    If the directory is invalid, it will return false.
    """

    if not is_valid_regulation_list_folder(input_path):
        return False

    info_filename = os.path.join(input_path, "info.csv")
    data_filename = os.path.join(input_path, "regulated_taxids.csv")
    child_lut_filename = os.path.join(input_path, "children_of_regulated_taxids.csv")

    _import_regulation_list_info(info_filename) # Always load first.
    _import_regulation_taxid_data(data_filename)
    _import_child_to_regulated_taxid_relationship(child_lut_filename)

    return True
