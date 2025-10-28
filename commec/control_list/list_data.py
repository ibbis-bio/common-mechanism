
import logging
import pandas as pd

from .containers import RegulationList

logger = logging.getLogger(__name__)

# Module Storage globals:

# The information of a single list, key on the list acronym.
REGULATION_LISTS : dict[str, RegulationList] = {}
# Information for every taxid within a list, keyed on list acronym, and taxid.
REGULATED_TAXID_ANNOTATIONS : pd.DataFrame = pd.DataFrame({
    "category": pd.Series(dtype="str"),
    "name": pd.Series(dtype="str"),
    "notes": pd.Series(dtype="str"),
    "derived_from": pd.Series(dtype="str"),
    "preferred_taxonomy_name": pd.Series(dtype="str"),
    "other_taxonomy_name": pd.Series(dtype="str"),
    "tax_id": pd.Series(dtype="Int64"),
    "genbank_protein": pd.Series(dtype="str"),
    "uniprot": pd.Series(dtype="str"),
    "list_acronym": pd.Series(dtype="str"),
    "target": pd.Series(dtype="str"),
    "hazard_group": pd.Series(dtype="str")
    })

# Map of children taxids to the regulated taxid in the lists.
CHILD_TAXID_MAP = pd.DataFrame({
    "child_taxid": pd.Series(dtype="Int64"),
    "regulated_taxid": pd.Series(dtype="Int64")
    })

def clear(target : str | None = None) -> bool:
    """
    Removes the targeted list from the module state, or
    if no target is provided, clears the entired module state.
    returns whether operation was successful.
    """
    global REGULATION_LISTS
    global REGULATED_TAXID_ANNOTATIONS

    if target and target in REGULATION_LISTS:
        # implement targetted removal logic.
        del REGULATION_LISTS[target]
        REGULATED_TAXID_ANNOTATIONS = REGULATED_TAXID_ANNOTATIONS[REGULATED_TAXID_ANNOTATIONS["list_acronym"] != target]
        # Consider updating the REG_TAXID_LISTS too.
        return True

    if not target:
        REGULATION_LISTS = {}
        REGULATED_TAXID_ANNOTATIONS = pd.DataFrame({
    "category": pd.Series(dtype="str"),
    "name": pd.Series(dtype="str"),
    "notes": pd.Series(dtype="str"),
    "derived_from": pd.Series(dtype="str"),
    "preferred_taxonomy_name": pd.Series(dtype="str"),
    "other_taxonomy_name": pd.Series(dtype="str"),
    "tax_id": pd.Series(dtype="Int64"),
    "genbank_protein": pd.Series(dtype="str"),
    "uniprot": pd.Series(dtype="str"),
    "list_acronym": pd.Series(dtype="str"),
    "target": pd.Series(dtype="str"),
    "hazard_group": pd.Series(dtype="str")
    })
        return True

    return False

def add_regulated_taxid_data(input_data : pd.DataFrame):
    """
    Calls concatenate to append new data to the regulated taxid annotations list.
    """
    global REGULATED_TAXID_ANNOTATIONS
    expected_cols = set(REGULATED_TAXID_ANNOTATIONS.columns)
    essential_cols = {"list_acronym", "tax_id"}

    # Check for missing columns
    if not essential_cols.issubset(input_data.columns):
    #if not (essential_cols in input_cols):
        raise ValueError(f"Input data must contain columns {essential_cols}, "
                         f"got {set(input_data.columns)}")

    input_data = input_data.reindex(columns=expected_cols)
    # Coerce invalid entries ("TBD", empty, etc.) to NaN, then cast to Int64
    input_data["tax_id"] = (
        pd.to_numeric(input_data["tax_id"], errors="coerce")
        .astype("Int64")
    )

    # The input data requires dtypes to be specified for any missing info.
    input_data = input_data.astype(REGULATED_TAXID_ANNOTATIONS.dtypes.to_dict())
    REGULATED_TAXID_ANNOTATIONS = pd.concat(
        [REGULATED_TAXID_ANNOTATIONS, input_data],
        ignore_index=True)

def add_child_lut_data(input_data: pd.DataFrame):
    """
    Append validated child TaxID mappings to CHILD_TAXID_MAP.
    Ensures correct columns, restricts to them, and removes duplicates.
    """
    global CHILD_TAXID_MAP
    expected_cols = set(CHILD_TAXID_MAP.columns)

    # Check for missing columns
    if not expected_cols.issubset(input_data.columns):
        raise ValueError(f"Input data must contain columns {expected_cols}, "
                         f"got {list(input_data.columns)}")

    # Restrict to only the expected columns
    input_data = input_data.reindex(columns=expected_cols)
    input_data["child_taxid"] = pd.to_numeric(input_data["child_taxid"], errors="coerce").astype("Int64")
    input_data["regulated_taxid"] = pd.to_numeric(input_data["regulated_taxid"], errors="coerce").astype("Int64")
    input_data = input_data.astype(CHILD_TAXID_MAP.dtypes.to_dict())
    CHILD_TAXID_MAP = pd.concat(
        [CHILD_TAXID_MAP, input_data], ignore_index=True
        ).drop_duplicates().reset_index(drop=True)

def add_regulated_list(new_list : RegulationList) -> bool:
    """
    Wrapper for safely adding a list to global data.
    """
    global REGULATION_LISTS

    # Overwrite protection:
    existing = REGULATION_LISTS.get(new_list.acronym)
    # Save the list.
    if existing:
        if existing == new_list:
            REGULATION_LISTS[new_list.acronym] = new_list
            return True
        logger.warning("Attempting to assign different list information to the same"
        "list identifier:\nExisting:\n%sNew:\n%s. List will be skipped.", existing, new_list)
        return False

    REGULATION_LISTS[new_list.acronym] = new_list
    return True
