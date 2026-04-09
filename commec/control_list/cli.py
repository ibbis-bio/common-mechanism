"""
Command-line-interface (cli) functionality for the Regulation module.
Argument declarations and 
"""
import os
import argparse
import pandas as pd
import importlib
from pathlib import Path
from commec.utils.file_utils import directory_arg, file_arg
import commec.config.yaml_io as YamlIO
from commec.config.constants import DEFAULT_CONFIG_YAML_PATH
from .containers import (
    ControlListOutput,
    ControlListContext,
)
from . import list_data as data

def add_args(parser_obj: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """
    Add Control List module arguments to an ArgumentParser object.
    """
    input_options = parser_obj.add_mutually_exclusive_group()
    input_options.add_argument(
        "-d",
        "--databases",
        dest="database_dir",
        type=directory_arg,
        #required = False,
        help="Path to parent directory containing Control List databases,"
        " the head of which should contain a region_definitions.json file",
    )
    input_options.add_argument(
        "-y",
        "--config",
        dest="yaml_file",
        type=file_arg,
        #required = True,
        help="Path to config file used by commec screen, to use the control list"
        " configuration as determined from that config yaml file.",
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
        action="store_true",
        help="Print a summary of all imported Control Lists",
    )
    parser_obj.add_argument(
        "-a",
        "--accessions",
        dest="showtaxids",
        nargs="+",
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

def format_control_lists(verbosity = False):
    """
    Summarises all loaded regulation list information,
    as well as their compliance under regional context.
    """
    output = ""
    if verbosity:
        output = "The following Control Lists apply: "
        for _, value in data.CONTROL_LISTS.items():
            number_of_regulated_taxids = (data.CONTROL_LIST_ANNOTATIONS["list_acronym"] == value.acronym).sum()
            output += f"\n{value.description()}\nRegulated Taxid Entries: {number_of_regulated_taxids}, Status : {value.status}"
        output += f"\n\n[Total number of Taxid Relationships:{data.ACCESSION_MAP.shape[0]}]"
        return output + "\n"

    # Table based output for reduced verbosity.
    rows = []
    for _, value in data.CONTROL_LISTS.items():
        rows.append({
            "Control List": value.name,# if len(value.name) < 50 else value.name[:50]+"...",
            "Acronym": value.acronym,
            "Region": value.region.acronym,
            "# Entries": (data.CONTROL_LIST_ANNOTATIONS["list_acronym"] == value.acronym).sum(),
            "Status": value.status
        })
    output = pd.DataFrame(rows, columns=["Control List", "Acronym", "# Entries","Region", "Status"]).to_string(index = False)
    return output

def format_control_list_annotation(input_data : list[ControlListOutput], input_context : list[ControlListContext]):
    """
    Returns a formatted string of the controlled taxid information for logging purposes.
    """
    plural = (len(input_data) > 1)
    output = "Controlled by the following lists:\n" if plural else ""
    offset = "       > " if plural else ""
    for i, output_info in enumerate(input_data):
        derived_string = input_context[i].derived_from or ""
        if derived_string:
            derived_string = derived_string + " "
        output += (offset + output_info.category + " "
                    + output_info.name
                    + f" controlled {derived_string}by {output_info.list}" + "\n")
    return output

def generate_output_summary_csv(output_filepath : str | os.PathLike):
    """
    Generates an output csv of the current controlled Annotations data
    imported into commec.
    Doesn't include invalid imported data with no Accession method.
    """

    output_data = data.CONTROL_LIST_ANNOTATIONS.copy(deep = True)
    
    # Currently removed as display_name, and list_item are new ways of working...
    #output_data["name"] = (
    #    output_data["preferred_taxonomy_name"].replace("", pd.NA)
    #    .combine_first(output_data["display_name"])
    #    .combine_first(output_data["other_taxonomy_name"].replace("", pd.NA))
    #)
    # output_data.drop(columns = ["preferred_taxonomy_name"], inplace=True)
    
    # Step 1: create dummy columns for list_acronym
    dummies = pd.get_dummies(output_data["list_acronym"], dtype = int)

    # Step 2: aggregate by index (since there may be multiple rows per index)
    indicators = dummies.groupby(output_data.index, sort=False).max()

    # Step 3: merge back with the *deduplicated* original dataframe
    # (dropping 'list_acronym' because it's now encoded)
    base = output_data.drop(columns="list_acronym").groupby(output_data.index, sort=False).first()
    result = base.join(indicators, how="outer")

    # Export - ensure .csv suffix
    output_path = Path(output_filepath)
    if output_path.suffix != ".csv":
        output_path = output_path.with_suffix(".csv")
    result.to_csv(output_path)

def read_config_yaml_for_control_list_info(config_yaml_filepath : os.PathLike | str):
    """
    Reads a config yaml, updated from the defaults, and parses the output for
    the control_list directory as per a commec screen run. Used instead of passing
    the directory of the control list
    
    :param config_yaml_filepath: Description
    :type config_yaml_filepath: os.PathLike | str
    """

    output_config = None

    # Read package-level configuration defaults
    default_yaml = importlib.resources.files("commec").joinpath(DEFAULT_CONFIG_YAML_PATH)
    if default_yaml.exists():
        output_config = YamlIO.load_config_from_yaml(str(default_yaml))
    else:
        raise FileNotFoundError(
            f"No default yaml found. Expected at {DEFAULT_CONFIG_YAML_PATH}"
            )

    # Override configuration with any in user-provided YAML file
    if os.path.exists(config_yaml_filepath):
        output_config = YamlIO.update_config_from_yaml(output_config, config_yaml_filepath)

    output_config = YamlIO.format_config_paths(output_config)

    return output_config
