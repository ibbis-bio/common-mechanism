# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

"""
Entry point for Control List, handles major API calls.

get_control_lists - Returns a list of the control lists information.
get_regulation - Returns control list information for an accession.
get_regulations - Batch version of get_regulation.
is_regulated - Returns truthiness for whether accession is regulated.
are_regulated - Batch version of is_regulated.
"""

import os
import logging
import argparse
from commec.utils.logger import setup_console_logging
from .containers import (
    Accession,
    derive_accession_format,
    ControlListOutput,
)
from . import list_data as __data
from . import initialisation as __init
from . import region as __region
from .cli import (
    add_args,
    format_control_lists,
    format_control_list_annotation,
    generate_output_summary_csv,
    read_config_yaml_for_control_list_info,
)

DESCRIPTION = """Tool for displaying information on
annotated control lists used during commec screen"""

logger = logging.getLogger(__name__)


def import_data(import_path: str | os.PathLike, regional_context: list[str] = None):
    """
    Entry point to load Control List data.
    Loads region definitions, then recursively loads data.
    """
    if not os.path.isdir(import_path):
        logger.warning("Provided path for control lists is invalid %s", import_path)
        return False

    # This needs to occur before we interpret regional context.
    __region.load_region_list_data(os.path.join(import_path, "region_definitions.json"))
    # Global, not control list specific ignore lists.
    __init.import_ignored_accesions(os.path.join(import_path, "ignored_taxids.csv"))

    # If the user puts in EU for example, we need to expand the set.
    cleaned_context = __region.get_regions_set(regional_context)
    logger.debug("Using the following regional context: \n")
    logger.debug(cleaned_context)

    # Recursively load the actual data.
    _import_data(import_path, cleaned_context)
    __init.tidy_control_list_data()

    if __data.CONTROL_LIST_ANNOTATIONS.size == 0:
        logger.warning("No control list annotations were found.")
        return False

    return True


def _import_data(import_path: str | os.PathLike, regional_context: list[str]):
    """
    Load a filepath recursively. Searches for valid "list folders"
    One found, the list folder is loaded into the modules state, and will
    """
    logger.debug("Checking path for list annotations: %s", import_path)
    if not os.path.isdir(import_path):
        logger.warning("Provided path invalid %s", import_path)
        return

    if __init.import_control_lists(import_path, regional_context):
        return

    logger.debug(
        "Invalid path: %s ... searching for more sub-directories...", import_path
    )

    # check for existance of sub folders and recurse on any present.
    for entry in os.scandir(import_path):
        if os.path.isdir(entry):
            _import_data(entry, regional_context)
    return


"""
Developers note:
It is tempting to mix should_ignore, and is_regulated, as they share some taxid
splitting logic and occur one after the other.
"""


def should_ignore(accession: str) -> bool:
    """
    Check whether this accession should simple be ignored by commec screen.
    """
    return str(accession) in __data.IGNORED_ACCESSION


def is_regulated(accession: str) -> bool:
    """
    Same as get_regulation, but optimised for speed — returns True/False
    for whether there is any control list data for the given accession.
    """

    accession_hash = Accession(accession)
    index_values = __data.CONTROL_LIST_ANNOTATIONS.index

    # Early exit if not present in control lists.
    in_map = accession_hash.code in __data.ACCESSION_MAP["child_taxid"].values
    in_index = accession_hash in index_values
    if not in_map and not in_index:
        return False
    return True


def get_cluster_hash(taxid: str) -> str | None:
    """
    Returns the highest parent, given the input accession.
    The highest parent is the child->controlled_taxid mapping that no longer
    maps to another taxid, whilst also being in the control list annotation index.
    Returns None if no regulated top-level ancestor exists.
    """

    # Early exit if not present in control lists.
    if not is_regulated(taxid):
        return None

    # Figure out the most parenty parent.
    index_values = __data.CONTROL_LIST_ANNOTATIONS.index
    accession_hash = Accession(taxid)
    accessions_to_check = {accession_hash}

    promote_child_to_parents = True
    while promote_child_to_parents:
        new_accession_set = set()
        promote_child_to_parents = False
        for accession in accessions_to_check:
            # Collect parent TaxIDs, if any
            taxid_parents = __data.ACCESSION_MAP.loc[
                __data.ACCESSION_MAP["child_taxid"] == accession.code,
                "controlled_taxid",
            ]
            if not taxid_parents.empty:
                new_accession_set.update(Accession(t) for t in taxid_parents)
                promote_child_to_parents = True
                continue
            # No further parent — this is a top-level node, keep it
            new_accession_set.add(accession)

        if accessions_to_check == new_accession_set:
            break
        accessions_to_check = new_accession_set
        promote_child_to_parents = True
        # if len(accessions_to_check) == 1:
        #    break

    if len(accessions_to_check) == 0:
        raise RuntimeError(
            f"get_cluster_hash({taxid!r}): traversal produced an empty set — "
            "taxid passed the early-exit check but all parents were consumed without resolution."
        )
    if len(accessions_to_check) > 1:
        logger.warning(
            "More than 1 highest parent for %s: %s", taxid, accessions_to_check
        )

    accession_to_check = next(iter(accessions_to_check))

    if accession_to_check not in index_values:
        logger.warning(
            "Accession parent (derived from %s) is apparently not in Control Lists %s",
            taxid,
            accession_to_check,
        )
        parent_accession = __data.ACCESSION_MAP.loc[
            __data.ACCESSION_MAP["child_taxid"] == accession_to_check.code,
            "controlled_taxid",
        ]
        logger.warning("Derived a parent: %s", parent_accession)

        return None

    return accession_to_check


def get_regulation(accession: str) -> list[ControlListOutput]:
    """
    Check the given Accession against all imported control lists.
    The input Accession can be a TaxID, GenBank protein, or Uniprot ID.
    If the input was a TaxID, any parent TaxIDS will also be recursively checked across
    all control lists.
    The output is a list of every regulation
    attributed to the original accession, in the form of
    a tuple containing the list info, as well as the
    taxid specific regulation information.
    """
    output_data: list[ControlListOutput] = []

    # Modify based on input accession format:
    accession_hash = Accession(accession)
    accession_to_check = [accession_hash]
    logger.debug("Fetching parents of '%s' [%s]", accession, type(accession))
    taxid_parents_to_check = __data.ACCESSION_MAP[
        __data.ACCESSION_MAP["child_taxid"] == accession_hash.code
    ]["controlled_taxid"].to_list()
    logger.debug("Found taxid parents: %s", str(taxid_parents_to_check))
    taxid_parents_to_check = [Accession(taxid) for taxid in taxid_parents_to_check]
    accession_to_check.extend(taxid_parents_to_check)

    logger.debug("Accessions to check: %s", accession_to_check)

    # Get Accessions of controlled interest:
    filtered_regulated_taxid_annotations = __data.CONTROL_LIST_ANNOTATIONS[
        __data.CONTROL_LIST_ANNOTATIONS.index.isin(accession_to_check)
    ]

    logger.debug(
        "Filtered Output DBS: %s", filtered_regulated_taxid_annotations.to_string()
    )

    # For each annotation, process its output, and context
    for hash_taxid, row in filtered_regulated_taxid_annotations.iterrows():
        is_child = accession_hash != hash_taxid
        child_name = ""
        if is_child:  # We want the childs name in this instance.
            child_name = __data.ACCESSION_MAP[
                __data.ACCESSION_MAP["child_taxid"] == str(accession_hash)
            ]["child_name"].iloc[0]

        output_data.append(
            ControlListOutput(
                row["display_name"],
                row["category"],
                row["list_acronym"],
                row["species"],
                row["genus"],
                row["list_item"],
                is_child,
                child_name,
            )
        )

    # Return the list pairs of output data and contexts.
    if len(output_data) > 0:
        logger.debug(
            "Checking %s [%s] for regulation resulted in %i annotations",
            accession_hash.get_format(),
            accession,
            len(output_data),
        )
    return output_data


def get_control_lists(list_acronym=None):
    """
    Simple retrieval for the 'list of Control lists' information.
    Optionally, pass a list acroynm, and retreive the information of that
    specific list.

    Returns None, or the ControlList, if list_acroynm was provided.
    Returns the list of all ControlLists if no input provided.
    """
    if list_acronym:
        return __data.CONTROL_LISTS.get(list_acronym)

    output = list(__data.CONTROL_LISTS.values())
    return output


def run(arguments: argparse.Namespace):
    """
    Run CLI for list printing etc.
    Currently also functions as a convenient test bed
    for ensuring the import scripts have worked.
    Tidy in future.
    """

    # Start logging to console
    log_level = logging.DEBUG if arguments.verbose else logging.INFO
    setup_console_logging(log_level)
    logger.info(
        " The Common Mechanism : List", extra={"no_prefix": True, "box_down": True}
    )

    logger.debug("Parsing input parameters... %s", arguments.database_dir)

    regions = arguments.regions or None

    if not (arguments.showlists or arguments.showtaxids or arguments.output_prefix):
        logger.error(
            "commec list requires --lists/-l, --accessions/-a, or --output_prefix/-o as input."
        )
        return 1

    database_location = None
    if arguments.database_dir:
        database_location = arguments.database_dir
    elif arguments.yaml_file:
        config = read_config_yaml_for_control_list_info(arguments.yaml_file)
        try:
            database_location = config["databases"]["control_lists"]["path"]
        except KeyError as e:
            logger.error(
                "Provided yaml input contained invaid control list path information. %s",
                e,
            )

    if not database_location:
        logger.error(
            "Provide the location of the control list database directory (-d)"
            " or location of yaml configuration file (-y) for commec list to import."
        )
        return 2

    import_data(database_location, regions)

    if arguments.showlists:
        logger.info(" *----------* CONTROL LISTS *----------* ")
        logger.info(format_control_lists(True), extra={"no_prefix": True, "cap": True})

    if arguments.showtaxids:
        logger.info(" *----------* CONTROLLED TAXIDS *----------* ")
        logger.info(
            "Regulation Annotations for supplied taxids (#%i):\n",
            len(arguments.showtaxids),
        )
        for accession in arguments.showtaxids:
            accession_format = derive_accession_format(accession)
            if not accession_format:
                logger.error(
                    "Could not determine the accession format for input: %s", accession
                )
                continue

            if should_ignore(accession):
                logger.info("   > Taxid %s is on the ignored list!", accession)
                continue

            outcome = get_regulation(accession)
            logger.info(
                "   > Taxid %s: %s",
                accession,
                (
                    format_control_list_annotation(outcome)
                    if len(outcome) > 0
                    else "No control list annotation found."
                ),
            )

    if arguments.output_prefix:
        logger.info(
            'Writing output list summary to "%s.csv" ...', arguments.output_prefix
        )
        generate_output_summary_csv(arguments.output_prefix)

    logger.info("", extra={"no_prefix": True, "box_up": True})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    args = parser.parse_args()
    run(args)
