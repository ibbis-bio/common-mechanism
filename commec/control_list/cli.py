"""
Command-line-interface (cli) functionality for the Regulation module.
Argument declarations and 
"""
import os
import argparse
import pandas as pd
from commec.utils.file_utils import directory_arg
from .containers import (
    ControlListInfo,
    ControlList,
)
from . import list_data as data

def add_args(parser_obj: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """
    Add Control List module arguments to an ArgumentParser object.
    """
    parser_obj.add_argument(
        "-d",
        "--databases",
        dest="database_dir",
        type=directory_arg,
        default=None,
        help="Path to parent directory containing Control List databases,"
        " the head of which should contain a region_definitions.json file",
    )
    parser_obj.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        help="Output additional debug logs.",
        default=False,
        action="store_true"
    )
    parser_obj.add_argument(
        "-l",
        "--list",
        dest="showlists",
        default=False,
        action="store_true",
        help="Print a summary of all imported Control Lists",
    )
    parser_obj.add_argument(
        "-a",
        "--accessions",
        dest="showtaxids",
        nargs="+",
        default=[],
        help="Display any available list information for the supplied Accessions.",
    )
    parser_obj.add_argument(
        "-r",
        "--regions",
        dest="regions",
        nargs="+",
        default=[],
        help="A list of countries or regions to add context to control list compliance",
    )
    parser_obj.add_argument(
        "-o",
        "--output_prefix",
        dest="output_prefix",
        help="Save a summary of all ingested information to an output fileset, with provided prefix",
        default="",
    )
    # --pretty?
    # --markdown? csv tsv etc

    return parser_obj

def format_control_lists():
    """
    Summarises all loaded regulation list information,
    as well as their compliance under regional context.
    """
    output = "The following Control Lists apply: "
    for _, value in data.CONTROL_LISTS.items():
        number_of_regulated_taxids = (data.CONTROL_LIST_ANNOTATIONS["list_acronym"] == value.acronym).sum()
        output += f"\n{value}\nRegulated Taxid Entries: {number_of_regulated_taxids}, Status : {value.status}"
    output += f"\n    [Total number of Taxid Relationships:{data.ACCESSION_MAP.shape[0]}]"
    return output + "\n"

def format_control_list_annotation(input_data : list[tuple[ControlList, ControlListInfo]]):
    """
    Returns a formatted string of the regulated taxid information for logging purposes.
    """
    plural = (len(input_data) > 1)
    output = "Regulated by the following lists:\n" if plural else ""
    offset = "   > " if plural else ""
    for reglist, annotations in input_data:
        output += (offset + annotations.category + " "
                    + annotations.name
                    + " regulated by " + reglist.name + f" [{reglist.acronym}]" + "\n")
    return output

def generate_output_summary_csv(output_filepath : str | os.PathLike):
    """
    Generates an output csv of the current Regulated Annotations data
    imported into commec.
    Doesn't include invalid imported data with no Accession method.
    """

    output_data = data.CONTROL_LIST_ANNOTATIONS.copy(deep = True)
    output_data["name"] = (
        output_data["preferred_taxonomy_name"]
        .combine_first(output_data["name"])
        .combine_first(output_data["other_taxonomy_name"])
    )
    output_data.drop(columns = ["preferred_taxonomy_name"], inplace=True)
    # Step 1: create dummy columns for list_acronym
    dummies = pd.get_dummies(output_data["list_acronym"], dtype = int)

    # Step 2: aggregate by index (since there may be multiple rows per index)
    indicators = dummies.groupby(output_data.index, sort=False).max()

    # Step 3: merge back with the *deduplicated* original dataframe
    # (dropping 'list_acronym' because it's now encoded)
    base = output_data.drop(columns="list_acronym").groupby(output_data.index, sort=False).first()
    result = base.join(indicators, how="outer")

    # Export
    result.to_csv(output_filepath)


