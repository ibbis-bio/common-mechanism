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

# Regions:
# Regions are always represented as a size 2 tuple, (name, acronym)

@dataclass
class Region:
    """
    Name and Acronym for regional information as either 
    a country, or wider a organisation:
    e.g. European Union, EU
    e.g. New Zealand, NZ
    """
    name = ""
    acronym = ""

    def __str__(self):
        return self.acronym

    def __repr__(self):
        return self.name + " ["+self.acronym+"]"


@dataclass
class RegulationList:
    """
    Contains the name, acronym, url, and affected regions for a regulated list.
    Will be output
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
    # full_list_name : str = ""
    # list_url : str = ""
    # region : str = ""

    taxonomy_category : str = ""
    taxonomy_name : str = ""
    notes : str = ""
    preferred_taxonomy_name : str = ""
    list_acronym : str = ""
    target : str = ""
    hazard_group : str = ""
    in_reg_taxids : str = ""
    parent_taxid : int = 0

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
    region_str : list[str] = field(default_factory=list(str))
    regulation_level : RegulationLevel = field(default_factory=RegulationLevel) # At what level this list is regulated.
    # toxicity : str = "no data" - TBD


# Module Storage globals:
# The information of a single list, key on the list acronym.
REG_LISTS : dict[str, RegulationList] = {}
# Information for every taxid within a list, keyed on list acronym, and taxid.
REG_TAXID_LISTS : dict[str, dict[int, TaxidRegulation]] = {}
# To be decided - if each TaxidRegulation simply contains a parent taxid,
# and bridging entries are used for each intermediary, then we might
# not require a global taxid map LUT...
TAXID_MAP : dict[int, tuple[int, list[int]]] = {}

def clear(target : str | None = None) -> bool:
    """
    Removes the targeted list from the module state, or
    if no target is provided, clears the entired module state.
    returns whether operation was successful.
    """
    if target and target in REG_LISTS:
        # implement targetted rmeoval logic.
        del REG_LISTS[target]
        return True

    if not target:
        REG_LISTS = None
        return True

    return False