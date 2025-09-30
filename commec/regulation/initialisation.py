# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

"""
Scripts that control the initialisation of the regulations
module state and imports regulation lists.

Import is handled recursively.

A valid list folder has the following layout:
```
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
```

Import goes by the following workflow:

"""
import os
import logging
import pandas as pd

from commec.regulation.containers import (
    RegulationList,
    ListMode,
    Region
    )

import commec.regulation.containers as rc
from commec.regulation.region import get_regions_set

logger = logging.getLogger(__name__)

def import_regulations(input_path : str | os.PathLike, regions_of_interest : set[str] = None) -> bool:
    """
    Given a directory, checks that the directory is valid,
    If the directory is valid, it will load the regulation data.
    If the directory is invalid, it will return false.

    ### Inputs:
    *input_path* : **PathLike | str**, Input directory to be searched.
    *regions_of_interest* : **set[str]**, 
        A set of alpha_2 country codes representing regions of interest. 
        Regulation lists which share at least one of these codes will be 
        marked for full compliance.
    
    ### outputs:
        **Boolean**: Returns True if the correct list files were detected, and
            a Regulated List and annotatations were imported.
    """
    info_filename = os.path.join(input_path, "list_info.csv")
    data_filename = os.path.join(input_path, "regulated_taxids.csv")
    child_lut_filename = os.path.join(input_path, "children_of_regulated_taxids.csv")

    # Check required files.
    if not os.path.isdir(input_path): return False
    for file in [info_filename, data_filename, child_lut_filename]:
        if not os.path.isfile(file): return False

    logger.debug("Importing regulation list from %s", input_path)
    _import_regulation_list_info(info_filename) # Always load first.
    _import_regulation_taxid_data(data_filename)
    _import_child_to_regulated_taxid_relationship(child_lut_filename)
    update_regional_context(regions_of_interest)
    return True


def update_regional_context(
    regions_of_interest : set[str] = None,
    alternative_mode = ListMode.CONDITIONAL_NUM
    ):
    """
    Updates all loaded regulation lists based on a regional context.
    Regulation lists which affect the same regions as the context are
    marked as requiring full compliance.
    Regions that do not affect regions of interest are labelled
    with the alternative mode, by default conditional based on number.

    ### inputs:
    
    *regions_of_interest* : **set[str]**, 
        A set of alpha_2 country codes representing regions of interest. 
        Regulation lists which share at least one of these codes will be 
        marked for full compliance.
    *alternative_mode* : **ListMode**, 
        What mode the regulation lists will be marked in the absence of 
        affecting any region of interest. Defaults to conditional based on 
        number of other lists also hit.
    
    ### outputs:
        *None*
    """

    if not regions_of_interest:
        regions_of_interest = set(["all"])

    force_all_regions = ("all" in regions_of_interest)

    if force_all_regions:
        logger.debug("Regulation list compliance set to affect all regions.")

    for reg_list in rc.REGULATION_LISTS.values():
        # Skip if forcing all regions.
        if force_all_regions:
            reg_list.status = ListMode.COMPLIANCE
            continue

        # Update regulation list mode based on region.
        list_affected_regions = get_regions_set(reg_list.regions)

        common_regions = set(list_affected_regions) & set(regions_of_interest)

        if common_regions:
            logger.debug("%s contains shared regions with context: %s",
                         reg_list.name, str(common_regions))
            reg_list.status = ListMode.COMPLIANCE
        else:
            logger.debug("%s no regions of context: %s",
                         reg_list.name, str(list_affected_regions))
            reg_list.status = alternative_mode


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
            [Region(name = row["region_name"],
                    acronym = row["region_code"])],
            ListMode.COMPLIANCE)

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
    rc.add_child_lut_data(child_lut)


def post_process_regulation_data():
    """
    These dataframes are always access via taxid, genbank, or uniprot,
    so it is much faster to index them based off this for querying.

    We also perform some data cleanup here, as well as take the opportunity to
    report to the user any tidying issues - such as duplicates with differing
    metadata.
    """
    rc.CHILD_TAXID_MAP.drop_duplicates()

    # Reindex the data based on the accession type column.
    rc.REGULATED_TAXID_ANNOTATIONS["accession"] = rc.REGULATED_TAXID_ANNOTATIONS.apply(
        lambda row: rc.Accession(
            taxid=row["tax_id"],
            genbank=row["genbank_protein"],
            uniprot=row["uniprot"], row = row),
        axis=1
    )

    # Duplicate annotatons may occur, but we only truely care about differences in the list acronym.
    # i.e. where a taxid is regulated from multiple sources.
    # We will still warn the user on differently formated duplicates.
    
    # Drop duplicates before indexing, using strict and non-strict strategy.
    bad_duplicates = rc.REGULATED_TAXID_ANNOTATIONS.drop_duplicates()
    rc.REGULATED_TAXID_ANNOTATIONS.drop_duplicates(
        subset=["list_acronym", "accession"], inplace=True)

    # Index and drop bad entries.
    rc.REGULATED_TAXID_ANNOTATIONS = rc.REGULATED_TAXID_ANNOTATIONS[
        rc.REGULATED_TAXID_ANNOTATIONS["accession"] != rc.Accession(taxid=0)]
    rc.REGULATED_TAXID_ANNOTATIONS.set_index("accession", inplace=True, drop = True)
    if rc.Accession(taxid=0) in rc.REGULATED_TAXID_ANNOTATIONS.index:
        rc.REGULATED_TAXID_ANNOTATIONS.drop(rc.Accession(taxid=0), inplace=True)
    
    # Remove bad entries, index for comparison.
    bad_duplicates = bad_duplicates[bad_duplicates["accession"] != rc.Accession(taxid=0)]
    bad_duplicates.set_index("accession", inplace=True, drop = False)

    # Merge comparison
    diff = pd.merge(bad_duplicates,
                    rc.REGULATED_TAXID_ANNOTATIONS,
                    how="left", indicator=True).query(
                        '_merge == "left_only"').drop(
                            columns="_merge")
    if not diff.empty:
        logger.warning("The following imported regulated annotations"
                       " were duplicates with differing metadata:\n%s",
                       diff[["accession", "name","category","list_acronym"]].to_string(index = False))

    logger.debug("Loaded the following regulation list dataset: Top 20:\n%s",
                 rc.REGULATED_TAXID_ANNOTATIONS.head(20).to_string())
