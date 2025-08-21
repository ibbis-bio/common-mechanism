# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

"""


"""
import os
import logging
import argparse
from commec.utils.logger import (
    setup_console_logging,
    setup_file_logging,
    set_log_level,
)
from commec.regulation.containers import TaxidRegulation, RegulationList
import commec.regulation.containers as data
from commec.regulation.initialisation import import_regulations
from commec.utils.file_utils import directory_arg, file_arg

DESCRIPTION = """Tool for displaying information on 
annotated regulated lists used during commec screen"""

logger = logging.getLogger(__name__)

def load_regulation_data(import_path : str | os.PathLike):
    """
    Load a filepath recursively. Searches for valid "list folders"
    One found, the list folder is loaded into the modules state, and will
    """
    logger.debug("Checking path for list annotations: %s", import_path)

    if import_regulations(import_path):
        return

    logger.debug("Invalid path: %s ... searching for more sub-directories...", import_path)

    # check for existance of sub folders and recurse on any present.
    for entry in os.scandir(import_path):
        if os.path.isdir(entry):
            load_regulation_data(entry)
    return

def get_regulation(taxid : int) -> list[tuple[RegulationList, TaxidRegulation]]:
    """
    Check the given TaxID against all imported regulated lists.
    The input taxid is 
    Any parent TaxIDS will also be recursively checked across 
    all regulated lists. The output is a list of every regulation 
    attributed to the original taxid and its parents, in the form of 
    a tuple containing the list, as well as the 
    taxid specific regulation information.
    """
    output_data : list[tuple[RegulationList, TaxidRegulation]] = []

    taxids_to_check = [taxid]
    taxid_parents_to_check = data.CHILD_TAXID_MAP[data.CHILD_TAXID_MAP["TaxID"] == taxid]["ParentTaxID"].to_list()
    taxids_to_check.append(taxid_parents_to_check)

    filtered_regulated_taxid_annotations = data.REGULATED_TAXID_ANNOTATIONS[
        data.REGULATED_TAXID_ANNOTATIONS["taxid"] in taxids_to_check]

    for _, row in filtered_regulated_taxid_annotations.iterrows():
        taxid_regulation_info = TaxidRegulation(
            row["taxonomy_category"],
            row["taxonomy_name"],
            row["notes"],
            row["preferred_taxonomy_name"],
            row["Taxid"],
            row["list_acronym"],
            row["target"],
            row["hazard_group"]
        )
        list_data = data.REGULATION_LISTS[taxid_regulation_info.list_acronym]
        output_data.append((list_data, taxid_regulation_info))

    return output_data

def print_regulation_list_information():
    ...


### Exact CLI arguments to be decided.
def add_args(parser_obj: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """
    Add module arguments to an ArgumentParser object.
    """
    parser_obj.add_argument(
        "-d",
        "--databases",
        dest="database_dir",
        type=directory_arg,
        default=None,
        help="Path to directory containing reference databases (e.g. taxonomy, protein, HMM)",
    )
    parser_obj.add_argument(
        "-y",
        "--config",
        dest="config_yaml",
        help="Configuration for screen run in YAML format, including custom database paths",
        default="",
    )

    parser_obj.add_argument(
        "-l",
        "--list",
        dest="showlists",
        default=False,
        action="store_true",
        help="Display annotation list information",
    )

    parser_obj.add_argument(
        "-t",
        "--taxids",
        dest="showtaxiddata",
        default=False,
        action="store_true",
        help="Display summary statistics on taxids",
    )



    # --pretty?

    # --markdown? 

    # --regions?

    return parser_obj


def run(args: argparse.Namespace):
    """Run CLI with an parsed argument parser input."""

    # Start logging to console
    log_level = logging.DEBUG
    setup_console_logging(log_level)
    logger.info(" The Common Mechanism : List", extra={"no_prefix": True, "box_down" : True})

    logger.debug("Parsing input parameters... %s", args.database_dir)

    if args.database_dir:
        logger.debug("Starting to load!")
        load_regulation_data(args.database_dir)

    logger.info("The following Regulation Lists have been identified: ")
    for _, value in data.REGULATION_LISTS.items():
        number_of_regulated_taxids = (data.REGULATED_TAXID_ANNOTATIONS["list_acronym"] == value.acronym).sum()
        logger.info("%s\nRegulated Taxid Entries: %s",value, number_of_regulated_taxids)

    logger.info("\nTotal number of Taxid Relationships: %i", data.CHILD_TAXID_MAP.shape[0])

    logger.debug("", extra={"no_prefix": True, "box_up" : True})

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    args = parser.parse_args()
    run(args)
