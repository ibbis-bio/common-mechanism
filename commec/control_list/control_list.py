# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

"""
Entry point for Control List, handles major API calls.

get_control_lists
get_regulation
is_regulated


"""
import os
import logging
import argparse
from commec.utils.logger import setup_console_logging
from .containers import (
    Accession,
    derive_accession_format,
    ControlListOutput,
    ControlListContext,
)
from . import list_data as __data
from . import initialisation as __init
from . import region as __region
from .cli import (
    add_args,
    format_control_lists,
    format_control_list_annotation,
    generate_output_summary_csv,
)

DESCRIPTION = """Tool for displaying information on
annotated regulated lists used during commec screen"""

logger = logging.getLogger(__name__)

def import_data(import_path : str | os.PathLike,
                         regional_context : list[str] = None):
    """
    Entry point to load Control List data.
    Loads region definitions, then recursively loads data.
    """
    # This needs to occur before we interpret regional context.
    __region.load_region_list_data(os.path.join(import_path, "region_definitions.json"))

    # If the user puts in EU for example, we need to expand the set.
    cleaned_context = __region.get_regions_set(regional_context)
    logger.debug("Using the following regional context: \n")
    logger.debug(cleaned_context)

    # Recursively load the actual data.
    _import_data(import_path, cleaned_context)
    __init.tidy_control_list_data()

def _import_data(import_path : str | os.PathLike,
                          regional_context : list[str]):
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

    logger.debug("Invalid path: %s ... searching for more sub-directories...", import_path)

    # check for existance of sub folders and recurse on any present.
    for entry in os.scandir(import_path):
        if os.path.isdir(entry):
            _import_data(entry, regional_context)
    return

def should_ignore(accession : str) -> bool:
    """
    Check whether this accession should simple be ignored by commec screen.
    """
    return (accession in __data.IGNORED_ACCESSION.to_numpy())

def is_regulated(accession : str) -> bool:
    """
    Same as get_regulation, but optimised for speed — returns True/False
    for whether there is any control list data for the given accession.
    """

    accession_hash = Accession(accession)
    accession_to_check = {accession_hash}

    # Collect parent TaxIDs, if any
    taxid_parents = __data.ACCESSION_MAP.loc[
        __data.ACCESSION_MAP["child_taxid"] == accession, "regulated_taxid"
    ]
    accession_to_check.update(Accession(taxid) for taxid in taxid_parents)

    # Check for intersection between sets of indexes vs accessions to check.
    index_values = __data.CONTROL_LIST_ANNOTATIONS.index
    return not accession_to_check.isdisjoint(index_values)

def get_regulation(accession : str) -> tuple[list[ControlListOutput], list[ControlListContext]]:
    """
    Check the given Accession against all imported regulated lists.
    The input Accession can be a TaxID, GenBank protein, or Uniprot ID.
    If the input was a TaxID, any parent TaxIDS will also be recursively checked across 
    all regulated lists. 
    The output is a list of every regulation
    attributed to the original accession, in the form of 
    a tuple containing the list info, as well as the 
    taxid specific regulation information.
    """
    output_data : list[ControlListOutput] = []
    output_context : list[ControlListContext] = []

    # Modify based on input accession format:
    accession_hash = Accession(accession)
    accession_to_check = [accession_hash]
    taxid_parents_to_check = __data.ACCESSION_MAP[
        __data.ACCESSION_MAP["child_taxid"] == accession]["regulated_taxid"].to_list()
    taxid_parents_to_check = [Accession(taxid) for taxid in taxid_parents_to_check]
    accession_to_check.extend(taxid_parents_to_check)

    logger.debug("Accesions to check: %s", accession_to_check)

    # Get Accessions of regulated interest:
    filtered_regulated_taxid_annotations = __data.CONTROL_LIST_ANNOTATIONS[
        __data.CONTROL_LIST_ANNOTATIONS.index.isin(accession_to_check)]
    
    logger.debug("Filtered Output DBS: %s", filtered_regulated_taxid_annotations.to_string())

    for hash_taxid, row in filtered_regulated_taxid_annotations.iterrows():
        output_data.append(ControlListOutput(row["name"],
                                            row["category"],
                                            row["list_acronym"]))
        output_context.append(ControlListContext(str(row["derived_from"]),
                                                 (accession != hash_taxid)))

    if len(output_data) > 0:
        logger.debug("Checking %s [%s] for regulation resulted in %i annotations",
                     accession_hash.get_format(), accession, len(output_data))

    return output_data, output_context

def get_control_lists(list_acronym = None):
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
    logger.info(" The Common Mechanism : List", extra={"no_prefix": True, "box_down" : True})

    logger.debug("Parsing input parameters... %s", arguments.database_dir)

    regions = arguments.regions or None

    if not (arguments.showlists or arguments.showtaxids):
        logger.error("commec list requires --lists/-l or --accessions/-a as input.")
        return 1

    if arguments.database_dir:
        import_data(arguments.database_dir, regions)
        
    if arguments.showlists:
        logger.info(" *----------* CONTROL LISTS *----------* ")
        logger.info(format_control_lists(True), extra={"no_prefix": True, "cap" : True})

    if arguments.showtaxids:
        logger.info(" *----------* REGULATED TAXIDS *----------* ")
        logger.info("Regulation Annotations "
                    "for supplied taxids (#%i):\n", len(arguments.showtaxids))
        for accession in arguments.showtaxids:
            accession_format = derive_accession_format(accession)
            if not accession_format:
                logger.error("Could not determine the accession format for input: %s", accession)
                continue
            outcome, outcome_context = get_regulation(accession)
            logger.info("   > Taxid %s: %s",accession,
                        format_control_list_annotation(outcome, outcome_context))

    if arguments.output_prefix:
        logger.info("Writing output list summary to \"%s.csv\" ...", arguments.output_prefix)
        generate_output_summary_csv(arguments.output_prefix)

    logger.info("", extra={"no_prefix": True, "box_up" : True})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    args = parser.parse_args()
    run(args)

