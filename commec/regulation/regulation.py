# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

"""
Responsible for the ingestion of annotated data from region based regulated
lists. 

"""
import os
import re
import logging
import argparse
import pandas as pd
from commec.utils.logger import setup_console_logging
from commec.regulation.containers import TaxidRegulation, RegulationList
import commec.regulation.containers as data
from commec.regulation.initialisation import import_regulations
from commec.regulation.region import load_region_list_data, get_regions_set
from commec.regulation.cli import (
    add_args,
    regulation_list_information,
    regulation_taxid_information
)

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

    # Recursively load the actual data.
    _load_regulation_data(import_path, cleaned_context)

    # These dataframes are always access via taxid, genbank, or uniprot,
    # so it is much faster to index them based off this for querying.
    data.REGULATED_TAXID_ANNOTATIONS["accession"] = data.REGULATED_TAXID_ANNOTATIONS.apply(
        lambda row: data.Accession(
            taxid=row["tax_id"],
            genbank=row["genbank_protein"],
            uniprot=row["uniprot"]),
        axis=1
    )
    data.REGULATED_TAXID_ANNOTATIONS.set_index("accession", inplace=True)

    logger.debug("Loaded the following dataset: Top 100:\n%s", data.REGULATED_TAXID_ANNOTATIONS.head(100).to_string())

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


def derive_accession_type(identifier : str) -> data.AccessionFormat | None:
    """
    Uses regex to determine whether the input
    string is a TaxID, Uniprot, or Genbank accesion.
    Returns None if it cannot be determined.
    """
    patterns = {
        data.AccessionFormat.UNIPROT: re.compile(
            r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
        ),
        data.AccessionFormat.TAXID: re.compile(
            r"^[0-9]+$"
        ),
        data.AccessionFormat.GENBANK: re.compile(
            r"^(?:[A-Z][0-9]{5}|[A-Z]{2}[0-9]{6,8})(?:\.\d+)?$"
        ),
    }

    for fmt, pattern in patterns.items():
        if pattern.fullmatch(identifier):
            return fmt

    return None


def get_regulation(accession : str, accession_fmt : data.AccessionFormat) -> list[tuple[RegulationList, TaxidRegulation]]:
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
    output_data : list[tuple[RegulationList, TaxidRegulation]] = []

    # Modify based on input accession format:
    accession_to_check = []
    match(accession_fmt):
        case data.AccessionFormat.TAXID:
            accession_to_check = [data.Accession(taxid = accession)]
            taxid_parents_to_check = data.CHILD_TAXID_MAP[
                data.CHILD_TAXID_MAP["child_taxid"] == accession]["regulated_taxid"].to_list()
            taxid_parents_to_check = [data.Accession(taxid=tid) for tid in taxid_parents_to_check]
            accession_to_check.extend(taxid_parents_to_check)

        case data.AccessionFormat.GENBANK:
            accession_to_check = [data.Accession(genbank = accession)]

        case data.AccessionFormat.UNIPROT:
            accession_to_check = [data.Accession(uniprot = accession)]

    logger.debug("accesions to check: %s", accession_to_check)

    # Get Accessions of regulated interest:
    filtered_regulated_taxid_annotations = data.REGULATED_TAXID_ANNOTATIONS[
        data.REGULATED_TAXID_ANNOTATIONS.index.isin(accession_to_check)]
    
    logger.debug("Filtered Output DBS: %s", filtered_regulated_taxid_annotations.to_string())

    for hash_taxid, row in filtered_regulated_taxid_annotations.iterrows():
        taxid_regulation_info = TaxidRegulation.from_row(row, hash_taxid)
        list_data = data.REGULATION_LISTS[taxid_regulation_info.list_acronym]
        output_data.append((list_data, taxid_regulation_info))

    if len(output_data) > 0:
        logger.debug("Checking %s [%s] for regulation resulted in %i annotations", accession_fmt, accession, len(output_data))

    return output_data

def testing_environment(output_prefix):
    print(data.REGULATED_TAXID_ANNOTATIONS.columns)
    extra_cols = data.REGULATED_TAXID_ANNOTATIONS[["tax_id","name"]].drop_duplicates()
    result = pd.crosstab(data.REGULATED_TAXID_ANNOTATIONS["tax_id"],
                         data.REGULATED_TAXID_ANNOTATIONS["list_acronym"]).astype(bool).astype(int)
    
    result = pd.merge(result, extra_cols, on="tax_id", how="left")

    output_filename = output_prefix + ".csv"
    # Write result to CSV

    pathogens = pd.read_csv("../IGSC/Pathogens.csv")
    toxins = pd.read_csv("../IGSC/Toxins.csv")

    result = pd.merge(result, pathogens, on="tax_id", how="outer")

    result.to_csv(output_filename, index=False, encoding="utf-8")

    print(pathogens.head(10))
    print(toxins.head(10))

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
            outcome = get_regulation(int(taxid), data.AccessionFormat.TAXID)

            logger.info("Taxid %s: %s",taxid,regulation_taxid_information(outcome))        

    logger.info("", extra={"no_prefix": True, "box_up" : True})

    # Early exit if no outputs requested.
    if not arguments.output_prefix:
        return
    
    testing_environment(arguments.output_prefix)

    unique_accession = data.REGULATED_TAXID_ANNOTATIONS.index.unique()
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    args = parser.parse_args()
    run(args)

