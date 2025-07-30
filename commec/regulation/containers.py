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
    Container for information of a given taxid.
    """
    taxonomy_category : str = ""
    taxonomy_name : str = ""
    notes : str = ""
    preferred_taxonomy_name : str = ""
    # full_list_name : str = ""
    list_acronym : str = ""
    # list_url : str = ""
    target : str = ""
    hazard_group : str = ""
    # region : str = ""
    in_reg_taxids : str = ""

@dataclass
class RegulationLevel(StrEnum):
    """
    At what level a regulation is for i.e. "regulated at the species level."
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
REG_LISTS : dict[str, RegulationList] = {}
REG_TAXID_LISTS : dict[str, dict[int, TaxidRegulation]] = {}
TAXID_MAP : dict[int, tuple[int, list[int]]] = {}
