# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

"""
Scripts that control the initialisation of the regulations
module state and imports regulation lists.

Import is handled recursively.

A valid list folder has the following layout:

regulated-lists/  		# This is the level at which pass to yaml / cli
├── uk-coshh/			# Arbitrary filename, contains 1x list.  
├── austgroup/          # The valid list folder layout:
│   ├── regulated_taxids.csv	# Taxid -> Parent, List annotations.
│   ├── regulated_taxids_and_children.csv
|   |		“
|   |		taxID, Name, ParentTaxID, ParentName
|   |		”
│   ├── regions.txt		# List of regions this list affects.
│   │   └── “
|   |		full_list_name,list_acronym,list_url,region_name, region,region
|   |		Special list, SL, www.sl.com,European Union, NZ, UK, NL, US
|   |		”
│   └── ...
└── ...

"""
import os

def _import_child_to_parent_relationship(input_file : str | os.PathLike):
    """

    """
    ...

def _is_valid_regulation_list_folder(input_path : str | os.PathLike) -> bool:
    """
    Checks if the supplies folder is a valid regulation list.
    i.e. contains the regulated_taxids.csv, regions_info.csv
    """
    ...

def _import_regulation_list(input_path : str | os.PathLike) -> bool:
    """
    Imports a folder containing a valid regulation list.
    Returns false if folder is invalid.
    """
    ...

def import_regulations(input_path : str | os.PathLike)