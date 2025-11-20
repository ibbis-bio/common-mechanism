# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

"""
Scripts that control the initialisation of the regulations
module state and imports control lists.

Import is handled recursively across arbitrarily nested folders.
The top level folder must contain a region_definitions.json.
A valid list folder may have the following layout:
```
control-lists/  		     # This is the level at which pass to yaml / cli
├── region_definitions.json  # Additional regional mapping information, i.e. EU
├── uk-coshh/			     # Arbitrary filename, contains 1x list.  
├── austgroup/               # The valid list folder layout:
│   ├── regulated_taxids.csv                          # Taxid List annotations.
│   ├── children_of_regulated_taxids.csv
|   |		“
|   |		child_taxid, regulated_taxid
|   |		”
│   ├── lift_info.csv		               # List of regions this list affects.
│   │   └── “
|   |		full_list_name,list_acronym,list_url,region_name, region_acronym
|   |		Special list, SL, www.sl.com, European Union, EU
|   |		”
│   └── ...
└── ...
```
"""
import os
import logging
import pandas as pd

from .containers import (
    CategoryType,
    ControlList,
    ListMode,
    ListUseAcronym,
    Region,
    Accession
    )

from . import list_data as ld
from .region import get_regions_set

logger = logging.getLogger(__name__)

def import_control_lists(
    input_path : str | os.PathLike,
    regions_of_interest : set[str] = None
    ) -> bool:
    """
    Given a directory, checks that the directory is valid,
    If the directory is valid, it will load the control list data.
    If the directory is invalid, it will return false.

    ### Inputs:
    *input_path* : **PathLike | str**, Input directory to be searched.
    *regions_of_interest* : **set[str]**, 
        A set of alpha_2 country codes representing regions of interest. 
        Control lists which share at least one of these codes will be 
        marked for full compliance.
    
    ### outputs:
        **Boolean**: Returns True if the correct list files were detected, and
            a Controlled List and annotatations were imported.
    """
    info_filename = os.path.join(input_path, "list_info.csv")
    data_filename = os.path.join(input_path, "regulated_taxids.csv")
    child_lut_filename = os.path.join(input_path, "children_of_regulated_taxids.csv")
    ignored_filename = os.path.join(input_path, "ignored_accessions.csv")

    # Check required files.
    if not os.path.isdir(input_path):
        return False
    for file in [info_filename, data_filename, child_lut_filename]:
        if not os.path.isfile(file):
            return False

    logger.debug("Importing Control List from %s", input_path)
    _import_control_list_info(info_filename) # Always load first.
    _import_control_list_annotations(data_filename)
    _import_accession_mappings(child_lut_filename)
    _import_ignored_accessions(ignored_filename)
    update_regional_context(regions_of_interest)
    return True


def update_regional_context(
    regions_of_interest : set[str] = None,
    alternative_mode = ListMode.CONDITIONAL_NUM
    ):
    """
    Updates all loaded Control lists based on a regional context.
    Control lists which affect the same regions as the context are
    marked as requiring full compliance.
    Regions that do not affect regions of interest are labelled
    with the alternative mode, by default conditional based on number.

    ### inputs:
    
    *regions_of_interest* : **set[str]**, 
        A set of alpha_2 country codes representing regions of interest. 
        Control lists which share at least one of these codes will be 
        marked for full compliance.
    *alternative_mode* : **ListMode**, 
        What mode the Control lists will be marked in the absence of 
        affecting any region of interest. Defaults to conditional based on 
        number of other lists also hit.
    
    ### outputs:
        *None*
    """

    if not regions_of_interest:
        regions_of_interest = set(["all"])

    force_all_regions = ("all" in regions_of_interest)

    if force_all_regions:
        logger.debug("Control list compliance set to affect all regions.")

    for reg_list in ld.CONTROL_LISTS.values():
        # Skip if forcing all regions.
        if force_all_regions:
            reg_list.status = ListMode.COMPLIANCE
            continue

        # Update control list mode based on region.
        list_affected_regions = get_regions_set(reg_list.region)

        common_regions = set(list_affected_regions) & set(regions_of_interest)

        if common_regions:
            logger.debug("%s contains shared regions with context: %s",
                         reg_list.name, str(common_regions))
            reg_list.status = ListMode.COMPLIANCE
        else:
            logger.debug("%s no regions of context: %s",
                         reg_list.name, str(list_affected_regions))
            reg_list.status = alternative_mode


def _import_control_list_info(input_path : str | os.PathLike):
    """
    Imports the info.csv file from a control list folder.
    Ensures that each control list information is appropriate keyed by
    acronym into the REG_LISTS data container.
    Ensures that existing control lists are not overwritten, and 
    warns the user if overwritting unique data (i.e. acroynym clash) has occured.
    """
    list_info = pd.read_csv(input_path, sep=",", quotechar='"', dtype = str)
    for _, row in list_info.iterrows():
        logger.debug("Parsing list information: %s", row)

        try:
            listuseacronym = ListUseAcronym(row["use"].strip())
        except ValueError as e:
            logger.error("Invalid value used for control list \"use\" column from file %s, error: %s ",
                         input_path, e)
            return

        try:
            new_region = Region(name = row["region_name"].strip(),
                    acronym = row["region_code"].strip())
        except ValueError as e:
            logger.error("Invalid values used for region definitions from file %s, error: %s ",
                         input_path, e)
            return

        new_list = ControlList(
            row["list_name"].strip(),
            row["list_acronym"].strip(),
            row["list_url"].strip(),
            new_region,
            ListMode.COMPLIANCE,
            listuseacronym
            )

        if ld.add_control_list(new_list):
            continue

        # Error adding list, likely a bad duplicate, report to user.
        list_key = row["list_acronym"]
        existing_list = ld.CONTROL_LISTS.get(list_key)
        logger.error("Error when importing lists from %s, shared acronym"
                    " between existing list %s and %s, New list will be skipped.",
                    input_path, existing_list, new_list)

def _import_control_list_annotations(input_path : str | os.PathLike):
    """
    Imports annotated controlled accession information from the regulated_taxids.csv file
    within a control list provided to commec.
    Each controlled accession is checked to ensure the control list is valid before
    inclusion.
    Concatenates the control taxid info into the global dataframe.
    """
    taxid_info = pd.read_csv(input_path, dtype = str)
    # We detect multiple list acroynms in the format "ABC, DEF, GHI"
    # Result: cells become lists like ['ABC', 'DEF', 'GHI']
    taxid_info["list_acronym"] = (
        taxid_info["list_acronym"]
        .astype(str)  # Ensure strings
        .str.split(",")  # Split on commas
        .apply(lambda x: [s.strip() for s in x])  # Strip whitespace
    )

    taxid_info["tax_id"] = (
        taxid_info["tax_id"]
        .astype(str)  # Ensure strings
        .str.split(",")  # Split on commas
        .apply(lambda x: [s.strip() for s in x])  # Strip whitespace
    )

    # "Explode" the lists into separate rows
    taxid_info = taxid_info.explode("list_acronym", ignore_index=True)
    taxid_info = taxid_info.explode("tax_id", ignore_index=True)

    # Only include data whose list acronym exists.
    mask = taxid_info["list_acronym"].apply(
        lambda list_key: ld.CONTROL_LISTS.get(list_key) is not None)
    valid_list_taxid_info = taxid_info[mask]

    # Warn about dropped rows
    dropped = taxid_info[~mask]
    if not dropped.empty:
        logger.warning("The following list acronyms were not valid from %s", input_path)
        logger.warning(dropped[["name","tax_id", "list_acronym"]].to_string())

    # Append the new list data:
    ld.add_control_list_annotations(valid_list_taxid_info)


def _import_accession_mappings(input_path : str | os.PathLike):
    """
    Imports child to controlled accession look up data data from the 
    children_of_regulated_taxids.csv file within a control list provided to commec.
    Concatenates the child LUT info into the global dataframe.
    """
    child_lut = pd.read_csv(input_path, dtype = str)
    ld.add_child_lut_data(child_lut)


def _import_ignored_accessions(input_path : str | os.PathLike):
    """
    Imports child to controlled accession look up data data from the 
    ignored_taxids.csv file within a control list provided to commec.
    Concatenates the ignored info into the global dataframe.
    """
    if os.path.isfile(input_path):
        ignored_data = pd.read_csv(input_path, dtype = str)
        ld.add_ignored_accession_data(ignored_data)


def tidy_control_list_data():
    """
    These dataframes are always access via taxid
    so it is much faster to index them based off this for querying.

    We also perform data cleanup here, as well as take the opportunity to
    report to the user any tidying issues - such as duplicates with differing
    metadata. In the ideal case, commec provided control annotations should not
    log any errors here.
    """
    ld.ACCESSION_MAP.drop_duplicates()
    ld.IGNORED_ACCESSION.drop_duplicates()

    # Reindex the data based on the accession type column.
    # Currently taxid is used as the only accession format.
    # Future versions will index based on accession more generally here.
    ld.CONTROL_LIST_ANNOTATIONS["accession"] = ld.CONTROL_LIST_ANNOTATIONS.apply(
        lambda row: Accession(accession=row["tax_id"]),
        axis=1
    )

    # Standardize the Category field.
    def safe_category(value):
        try:
            return CategoryType.NONE if pd.isna(value) else CategoryType(value)
        except ValueError:
            return CategoryType.NONE

    ld.CONTROL_LIST_ANNOTATIONS["category"] = ld.CONTROL_LIST_ANNOTATIONS["category"].map(safe_category)

    # Report errors for bad entries
    bad_entries = ld.CONTROL_LIST_ANNOTATIONS[
        ld.CONTROL_LIST_ANNOTATIONS["accession"] == Accession(None)]
    bad_entries.loc[bad_entries["name"].str.len() >= 57, "name"] = (
        bad_entries["name"].str[:57].str.strip() + "..."
    )
    if not bad_entries.empty:
        logger.error("%i imported control list annotations"
                       " were bad entries with no TaxID Accession:\n%s"
                       "\n Run in --verbose mode for raw row input details.",
                       len(bad_entries.index),
                       bad_entries[["name","category","list_acronym"]].to_string(index = False))

    # Drop duplicates before indexing, using strict and non-strict strategy.
    bad_duplicates = ld.CONTROL_LIST_ANNOTATIONS.drop_duplicates()
    ld.CONTROL_LIST_ANNOTATIONS.drop_duplicates(
        subset=["list_acronym", "accession"], inplace=True)

    # Index and drop bad entries.
    ld.CONTROL_LIST_ANNOTATIONS = ld.CONTROL_LIST_ANNOTATIONS[
        ld.CONTROL_LIST_ANNOTATIONS["accession"] != Accession(None)]
    ld.CONTROL_LIST_ANNOTATIONS.set_index("accession", inplace=True, drop = True)
    if Accession(None) in ld.CONTROL_LIST_ANNOTATIONS.index:
        ld.CONTROL_LIST_ANNOTATIONS.drop(Accession(None), inplace=True)

    # Remove bad entries, index for comparison.
    bad_duplicates = bad_duplicates[bad_duplicates["accession"] != Accession(None)]
    bad_duplicates.set_index("accession", inplace=True, drop = False)

    # Merge comparison
    diff = pd.merge(bad_duplicates,
                    ld.CONTROL_LIST_ANNOTATIONS,
                    how="left", indicator=True).query(
                        '_merge == "left_only"').drop(
                            columns="_merge")
    if not diff.empty:
        logger.debug("The following imported control list annotations"
                       " were duplicates with differing metadata:\n%s",
                       diff[["accession","name","category","list_acronym"]].to_string(index = False))

    logger.debug("Loaded the following control list dataset: Top 20:\n%s",
                 ld.CONTROL_LIST_ANNOTATIONS.head(20).to_string())

