# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

"""


"""
import os
import logging
import argparse
from commec.utils.logger import (
    setup_console_logging
)
from commec.regulation.containers import TaxidRegulation, RegulationList
import commec.regulation.containers as data
from commec.regulation.initialisation import import_regulations
from commec.regulation.region import load_region_list_data, get_regions_set
from commec.utils.file_utils import directory_arg

DESCRIPTION = """Tool for displaying information on
annotated regulated lists used during commec screen"""

logger = logging.getLogger(__name__)

def load_regulation_data(import_path : str | os.PathLike,
                         regional_context : list[str] = None):
    """
    Entry point to load regulation data.
    Loads region definitions, then recursively loads data.
    """
    # This needs to occur before we interpret regional context.
    load_region_list_data(os.path.join(import_path, "region_definitions.json"))

    # If the user puts in EU for example, we need to expand the set.
    cleaned_context = get_regions_set(regional_context)
    logger.debug("Using the following regional context: \n")
    logger.debug(cleaned_context)

    # Load the actual data.
    _load_regulation_data(import_path, cleaned_context)

def _load_regulation_data(import_path : str | os.PathLike,
                          regional_context : list[str]):
    """
    Load a filepath recursively. Searches for valid "list folders"
    One found, the list folder is loaded into the modules state, and will
    """
    logger.debug("Checking path for list annotations: %s", import_path)

    if import_regulations(import_path, regional_context):
        return

    logger.debug("Invalid path: %s ... searching for more sub-directories...", import_path)

    # check for existance of sub folders and recurse on any present.
    for entry in os.scandir(import_path):
        if os.path.isdir(entry):
            _load_regulation_data(entry, regional_context)
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

    TODO: Update to accept Uniprot and genbank accessions, not just taxid.
    """
    #logger.debug("Checking taxid [%i] for regulation ", taxid)
    output_data : list[tuple[RegulationList, TaxidRegulation]] = []

    taxids_to_check = [taxid]
    taxid_parents_to_check = data.CHILD_TAXID_MAP[data.CHILD_TAXID_MAP["child_taxid"] == taxid]["regulated_taxid"].to_list()
    taxids_to_check.extend(taxid_parents_to_check)
    #logger.debug("Additional taxids to check: %s", taxid_parents_to_check)

    filtered_regulated_taxid_annotations = data.REGULATED_TAXID_ANNOTATIONS[
        data.REGULATED_TAXID_ANNOTATIONS["taxid"].isin(taxids_to_check)
    ]

    #logger.debug("Filtered Output DBS: %s", filtered_regulated_taxid_annotations.to_string())

    for _, row in filtered_regulated_taxid_annotations.iterrows():
        taxid_regulation_info = TaxidRegulation(
            row["taxonomy_category"],
            row["taxonomy_name"],
            row["notes"],
            row["preferred_taxonomy_name"],
            row["list_acronym"],
            int(row["taxid"]),
            row["target"],
            row["hazard_group"]
        )
        list_data = data.REGULATION_LISTS[taxid_regulation_info.list_acronym]
        output_data.append((list_data, taxid_regulation_info))

    if len(output_data) > 0:
        logger.debug("Checking taxid [%i] for regulation resulted in %i annotations", taxid, len(output_data))

    return output_data

def regulation_list_information():
    """
    Summarises all loaded regulation list information,
    as well as their compliance under regional context.
    """
    output = "The following Regulation Lists have been identified: "
    for _, value in data.REGULATION_LISTS.items():
        number_of_regulated_taxids = (data.REGULATED_TAXID_ANNOTATIONS["list_acronym"] == value.acronym).sum()
        output += f"\n{value}\nRegulated Taxid Entries: {number_of_regulated_taxids}, Status : {value.status}"
    output += f"\n    [Total number of Taxid Relationships:{data.CHILD_TAXID_MAP.shape[0]}]"
    return output + "\n"

def regulation_taxid_information(input_data : list[tuple[RegulationList, TaxidRegulation]]):
    """
    Returns a formatted string of the regulated taxid information for logging purposes.
    """
    plural = (len(input_data) > 1)
    output = "Regulated by the following lists:\n" if plural else ""
    offset = "   > " if plural else ""
    for reglist, annotations in input_data:
        output += (offset + annotations.taxonomy_category + " "
                    + annotations.preferred_taxonomy_name
                    + " regulated by " + reglist.name + f" [{reglist.acronym}]" + "\n")
    return output

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
        "-v",
        "--verbose",
        dest="verbose",
        help="Output debug logs.",
        default=False,
        action="store_true"
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
        dest="showtaxids",
        nargs="+",
        default=[],
        help="Display any available list information for the supplied taxids.",
    )
    parser_obj.add_argument(
        "-r",
        "--regions",
        dest="regions",
        nargs="+",
        default=[],
        help="A list of countries or regions to add context to list compliance",
    )

    # --pretty?
    # --markdown? csv tsv etc

    return parser_obj


def run(arguments: argparse.Namespace):
    """Run CLI for list printing etc. Currently also functions as a convenient test bed
    for ensuring the import scripts have worked. Tidy in future."""

    # Start logging to console
    log_level = logging.DEBUG if arguments.verbose else logging.INFO
    setup_console_logging(log_level)
    logger.info(" The Common Mechanism : List", extra={"no_prefix": True, "box_down" : True})

    logger.debug("Parsing input parameters... %s", arguments.database_dir)

    regions = arguments.regions or None

    if arguments.database_dir:
        load_regulation_data(arguments.database_dir, regions)

    if arguments.showlists:
        logger.info(" *----------* REGULATION LISTS *----------* ")
        logger.info(regulation_list_information())

    if arguments.showtaxids:
        logger.info(" *----------* REGULATED TAXIDS *----------* ")
        logger.info("Regulation Annotations "
                    "for supplied taxids (#%i):\n", len(arguments.showtaxids))
        for taxid in arguments.showtaxids:
            outcome = get_regulation(int(taxid))

            logger.info("Taxid %s: %s",taxid,regulation_taxid_information(outcome))

    logger.info("", extra={"no_prefix": True, "box_up" : True})

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    args = parser.parse_args()
    run(args)
