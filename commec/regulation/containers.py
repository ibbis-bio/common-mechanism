# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

"""
Defines Containers to show default layouts for regulation data.
What is the shape of the regulation data, and where to store it.
-----
    RegulationContainer : 
        A standard output that may be associated with a hit.
        Contains list, region, and level information.
    RegulationLevel :
        Determines valid level information strings.
    TaxidRegulationContainer :
        Storage of all regulation data for a given taxid from a given list.

-----
Also defines the storage data for 
    * regulation lists, 
    * per taxid regulation information per list, 
    * and child to parent taxid mapping.
"""
from dataclasses import dataclass, field
from enum import StrEnum
import pandas as pd

@dataclass
class Region:
    """
    Name and Acronym for regional information as either 
    a country, or wider a organisation:
    e.g. European Union, EU
    e.g. New Zealand, NZ
    """
    name : str = ""
    acronym : str = ""

    def __str__(self):
        return self.acronym

    def __repr__(self):
        return self.name + " ["+self.acronym+"]"


@dataclass
class RegulationList:
    """
    Contains the name, acronym, url, and affected regions for a regulated list.
    """
    name : str = ""
    acronym : str = ""
    url : str = ""
    regions : list[Region] = field(default_factory=list[Region])

@dataclass
class TaxidRegulation:
    """
    Container for regulatory information of a given taxid from
    a specific list.
    """
    taxonomy_category : str = ""
    taxonomy_name : str = ""
    notes : str = ""
    preferred_taxonomy_name : str = ""
    list_acronym : str = ""
    target : str = ""
    hazard_group : str = ""
    in_reg_taxids : str = ""

@dataclass
class RegulationLevel(StrEnum):
    """
    At what level a regulation is for i.e. "regulated at the species level."
    ??? These delineations require further thought before being practical.
    Specifically, protein may only apply to biorisk hits, not taxonomy screening.
    """
    ORGANISM = "organism" # All species of this organism are regulated.
    SPECIES = "species" # This specific species is regulated.
    PROTEIN = "protein" # This protein is regulated, regardless of organism.

@dataclass
class RegulationOutput:
    """
    Container for regulation list information. Formatted for use in output JSON.
    -----
    * `name` ( str ) : Common name identification for this regulation list.
    * `region` ( str ) : 2 Letter Country codes for which this regulation applies
    * `regulation_level` ( str ) : At what level does this regulation apply to the hit (e.g. organism, species, protein)
    -----
    See results.py for other examples.
    """
    name_str : str = ""
    region_str : list[str] = field(default_factory=list[str])
    # regulation_level : RegulationLevel = field(default_factory=RegulationLevel) # At what level this list is regulated.
    # toxicity : str = "no data" - TBD


# Module Storage globals:
# The information of a single list, key on the list acronym.
REGULATION_LISTS : dict[str, RegulationList] = {}
# Information for every taxid within a list, keyed on list acronym, and taxid.
# REG_TAXID_LISTS : dict[str, dict[int, TaxidRegulation]] = {}
REGULATED_TAXID_ANNOTATIONS : pd.DataFrame = pd.DataFrame(columns = [
    "taxonomy_category",
    "taxonomy_name",
    "notes",
    "preferred_taxonomy_name",
    "taxid",
    "list_acronym",
    "target",
    "hazard_group"
    ])
# Map of children taxids to the regulated taxid in the lists.
CHILD_TAXID_MAP : pd.DataFrame = pd.DataFrame(columns = [
    "TaxID",
    "ParentTaxID"
    ])

def clear(target : str | None = None) -> bool:
    """
    Removes the targeted list from the module state, or
    if no target is provided, clears the entired module state.
    returns whether operation was successful.
    """
    global REGULATION_LISTS

    if target and target in REGULATION_LISTS:
        # implement targetted removal logic.
        del REGULATION_LISTS[target]
        # Consider updating the REG_TAXID_LISTS too.
        return True

    if not target:
        REGULATION_LISTS = None
        return True

    return False

def add_regulated_taxid_data(input_data : pd.DataFrame):
    """
    Calls concatenate to append new data to the regulated taxid annotations list.
    """
    global REGULATED_TAXID_ANNOTATIONS
    expected_cols = {"taxonomy_category",
                    "taxonomy_name",
                    "notes",
                    "preferred_taxonomy_name",
                    "taxid",
                    "list_acronym",
                    "target",
                    "hazard_group"}

    # Check for missing columns
    if not expected_cols.issubset(input_data.columns):
        raise ValueError(f"Input data must contain columns {expected_cols}, "
                         f"got {list(input_data.columns)}")

    # Restrict to only the expected columns
    input_data = input_data[["taxonomy_category",
                    "taxonomy_name",
                    "notes",
                    "preferred_taxonomy_name",
                    "taxid",
                    "list_acronym",
                    "target",
                    "hazard_group"]]

    # Concatenate and drop exact duplicates
    REGULATED_TAXID_ANNOTATIONS = pd.concat(
        [REGULATED_TAXID_ANNOTATIONS, input_data], 
        ignore_index=True)

def add_child_lut_data(input_data: pd.DataFrame):
    """
    Append validated child TaxID mappings to CHILD_TAXID_MAP.
    Ensures correct columns, restricts to them, and removes duplicates.
    """
    global CHILD_TAXID_MAP

    expected_cols = {"TaxID", "ParentTaxID"}

    # Check for missing columns
    if not expected_cols.issubset(input_data.columns):
        raise ValueError(f"Input data must contain columns {expected_cols}, "
                         f"got {list(input_data.columns)}")

    # Restrict to only the expected columns
    input_data = input_data[["TaxID", "ParentTaxID"]]

    # Concatenate and drop exact duplicates
    CHILD_TAXID_MAP = (
        pd.concat([CHILD_TAXID_MAP, input_data], ignore_index=True)
          .drop_duplicates()
          .reset_index(drop=True)
    )

def add_regulated_list(input : RegulationList):
    global REGULATION_LISTS
    REGULATION_LISTS[input.acronym] = input
