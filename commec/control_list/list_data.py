"""
Storage for the module globals, containing both the
 * Control Lists
 * Control list annotations
 * Accession Map
 * Ignored Accessions
As well as methods for interacting safely with these globals:
 * clear()
 *
"""

import logging
import pandas as pd

from .containers import ControlList

logger = logging.getLogger(__name__)

# The information of the Control Lists, Indexed on the list acronym.
CONTROL_LISTS: dict[str, ControlList] = {}

# Information for every Accession within a list,
# unique for list acronym, and taxid.
# Indexed on Accession
CONTROL_LIST_ANNOTATIONS: pd.DataFrame = pd.DataFrame(
    {
        "list_acronym": pd.Series(dtype="str"),
        "list_item": pd.Series(dtype="str"),
        "tax_id": pd.Series(dtype="str"),
        "display_name": pd.Series(dtype="str"),
        "category": pd.Series(dtype="str"),
        "kingdom": pd.Series(dtype="str"),
        "genus": pd.Series(dtype="str"),
        "species": pd.Series(dtype="str"),
        "strain": pd.Series(dtype="str"),
    }
)

# Precalculated and imported map of child
# accessions to control list Accessions.
ACCESSION_MAP = pd.DataFrame(
    {
        "child_taxid": pd.Series(dtype="str"),
        "controlled_taxid": pd.Series(dtype="str"),
        "child_name": pd.Series(dtype="str"),
    }
)

# Precalculated and imported set of accessions
# Typically unclassified cousins of controlled accessions,
# And contains all vaccine and synthetic taxids.
IGNORED_ACCESSION: set[str] = set()


def clear(target: str | None = None) -> bool:
    """
    Removes the targeted list from the module state, or
    if no target is provided, clears the entired module state.
    returns whether operation was successful.
    """
    global CONTROL_LISTS
    global CONTROL_LIST_ANNOTATIONS
    global ACCESSION_MAP
    global IGNORED_ACCESSION

    if target and target in CONTROL_LISTS:
        # implement targetted removal logic.
        del CONTROL_LISTS[target]
        CONTROL_LIST_ANNOTATIONS = CONTROL_LIST_ANNOTATIONS[
            CONTROL_LIST_ANNOTATIONS["list_acronym"] != target
        ]
        return True

    if not target:
        CONTROL_LISTS = {}
        ACCESSION_MAP = pd.DataFrame(
            {
                "child_taxid": pd.Series(dtype="str"),
                "controlled_taxid": pd.Series(dtype="str"),
                "child_name": pd.Series(dtype="str"),
            }
        )
        IGNORED_ACCESSION = set()
        CONTROL_LIST_ANNOTATIONS = pd.DataFrame(
            {
                "list_acronym": pd.Series(dtype="str"),
                "list_item": pd.Series(dtype="str"),
                "tax_id": pd.Series(dtype="str"),
                "display_name": pd.Series(dtype="str"),
                "category": pd.Series(dtype="str"),
                "kingdom": pd.Series(dtype="str"),
                "genus": pd.Series(dtype="str"),
                "species": pd.Series(dtype="str"),
                "strain": pd.Series(dtype="str"),
            }
        )
        return True

    return False


def add_control_list(new_list: ControlList) -> bool:
    """
    Wrapper for safely adding a Control List to this modules global data.
    Handles List information collisions, returns True on success.
    """
    global CONTROL_LISTS

    # Overwrite protection:
    existing = CONTROL_LISTS.get(new_list.acronym)
    if existing:
        if existing == new_list:
            return True
        return False

    # Save the list as new:
    CONTROL_LISTS[new_list.acronym] = new_list
    return True


def add_control_list_annotations(input_data: pd.DataFrame):
    """
    Appends new annotation information, and controls for data collisions,
    only includes data as part of the expected columns.
    """
    global CONTROL_LIST_ANNOTATIONS
    expected_cols = set(CONTROL_LIST_ANNOTATIONS.columns)
    essential_cols = {"list_acronym", "tax_id"}

    # Check for missing columns
    if not essential_cols.issubset(input_data.columns):
        # if not (essential_cols in input_cols):
        raise ValueError(
            f"Input data must contain columns {essential_cols}, "
            f"got {set(input_data.columns)}"
        )

    input_data = input_data.reindex(columns=expected_cols)

    # The input data requires dtypes to be specified for any missing info.
    input_data = input_data.astype(CONTROL_LIST_ANNOTATIONS.dtypes.to_dict()).fillna("")

    CONTROL_LIST_ANNOTATIONS = pd.concat(
        [CONTROL_LIST_ANNOTATIONS, input_data], ignore_index=True
    )


def add_child_lut_data(input_data: pd.DataFrame):
    """
    Append pre-calculated child accession mappings to ACCESSION_MAP.
    Ensures correct columns, restricts to them, and removes duplicates.
    """
    global ACCESSION_MAP
    expected_cols = set(ACCESSION_MAP.columns)

    # Check for missing columns
    if not expected_cols.issubset(input_data.columns):
        raise ValueError(
            f"Input data must contain columns {expected_cols}, "
            f"got {list(input_data.columns)}"
        )

    # Restrict to only the expected columns
    input_data = input_data.reindex(columns=expected_cols)
    input_data = input_data.astype(ACCESSION_MAP.dtypes.to_dict())
    ACCESSION_MAP = (
        pd.concat([ACCESSION_MAP, input_data], ignore_index=True)
        .drop_duplicates()
        .reset_index(drop=True)
    )


def add_ignored_accession_data(input_data: pd.DataFrame):
    """
    Append pre-calculated ignored accessions to IGNORED_ACCESSION.
    """
    global IGNORED_ACCESSION

    if "ignored_taxid" not in input_data.columns:
        raise ValueError(
            f"Input data must contain column 'ignored_taxid', "
            f"got {list(input_data.columns)}"
        )

    new_accessions = set(input_data["ignored_taxid"].astype(str))
    IGNORED_ACCESSION.update(new_accessions)
    logger.debug(
        "Added %i ignored accessions. Total: %i",
        len(new_accessions),
        len(IGNORED_ACCESSION),
    )
