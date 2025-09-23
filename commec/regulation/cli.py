import argparse
from commec.utils.file_utils import directory_arg
from commec.regulation.containers import TaxidRegulation, RegulationList
import commec.regulation.containers as data


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
    parser_obj.add_argument(
        "-o",
        "--output_prefix",
        dest="output_prefix",
        help="Save all pertinant information to an output fileset, with provided prefix",
        default="",
    )
    # --pretty?
    # --markdown? csv tsv etc

    return parser_obj

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
        output += (offset + annotations.category + " "
                    + annotations.preferred_taxonomy_name
                    + " regulated by " + reglist.name + f" [{reglist.acronym}]" + "\n")
    return output