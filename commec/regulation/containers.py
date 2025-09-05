# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

"""
Defines Containers to show default layouts for regulation data.
What is the shape of the regulation data, and where to store it.

    RegulationContainer : 
        A standard output that may be associated with a hit.
        Contains list, region, and level information.
    RegulationLevel :
        Determines valid level information strings.
    TaxidRegulationContainer :
        Storage of all regulation data for a given taxid from a given list.

Also defines the storage data for 
    * regulation lists, 
    * per taxid regulation information per list, 
    * and child to parent taxid mapping.
"""
import logging
from dataclasses import dataclass, field, fields
from typing import Optional
from enum import StrEnum
import pandas as pd

logger = logging.getLogger(__name__)

class AccesionFormat(StrEnum):
    """
    Supported hashable accesion formats.
    """
    TAXID = "Taxonomy ID"
    GENBANK = "Genbank Protein Accession"
    UNIPROT = "Uniprot ID"

    def __repr__(self):
        return f"<Fmt:{self}>"

@dataclass(frozen=True)
class Accession:
    """
    A unified hashable identifier for 
    taxonomy IDs, GenBank, or UniProt accessions.

    Used for hashed index accession for rapid annotations retrieval
    from a pandas DataFrame.
    """
    code: str
    type: AccesionFormat   # one of: "taxid", "genbank", "uniprot"

    def __init__(self, taxid: Optional[int]=None,
                       uniprot: Optional[str]=None,
                       genbank: Optional[str]=None):
        
        # ensure only one is set
        values = [(taxid, AccesionFormat.TAXID),
                  (uniprot, AccesionFormat.UNIPROT),
                  (genbank, AccesionFormat.GENBANK)]
        non_empty = [(v, t) for v, t in values if not pd.isna(v)]
        if len(non_empty) != 1:
            #raise ValueError(f"Accession must have exactly one non-empty value, got {non_empty} from {values}")
            non_empty = [("0", AccesionFormat.TAXID)]


        object.__setattr__(self, "code", str(non_empty[0][0]))
        object.__setattr__(self, "type", non_empty[0][1])

    def get(self) -> tuple[str, str]:
        """Return (code, type)."""
        return self.code, self.type

    #def __hash__(self):
    #    # frozen dataclasses already hash on all fields, but explicit is nice
    #    return hash((self.code, self.type))

    def __repr__(self):
        return f"<Accession {self.type}:{self.code}>"

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
        return self.name

    def __repr__(self):
        return self.name + " ["+self.acronym+"]"

class ListMode(StrEnum):
    """
    Describes how a list is interpreted, based on regional context.
    When a Taxonomy ID is identified as part of a list, whether or not we treat
    that taxid as regulated or not will depend on the following mode.
    
    * COMPLIANCE - All taxids from this list are regulated.
    * CONDITIONAL_NUM - Only mark taxid as regulated if it appears in more than 1 list.
    * COMPLIANCE_WARN - Mark Taxid as only a WARNING.
    * IGNORE - Ignore this list entirely.
    """
    COMPLIANCE = "Compliance"
    CONDITIONAL_NUM = "Conditional Compliance"
    COMPLIANCE_WARN = "Comply with Warning"
    IGNORE = "Ignored"

@dataclass
class RegulationList:
    """
    Contains the name, acronym, url, and affected regions for a regulated list.
    """
    name : str = ""
    acronym : str = ""
    url : str = ""
    regions : list[Region] = field(default_factory=list[Region])
    status : ListMode = field(default_factory=ListMode)

    def __str__(self):
        regions_text = ""
        for r in self.regions:
            regions_text += str(r) + ", "
        regions_text = regions_text[:-2]
        return f"[{self.acronym}] {self.name} - {regions_text}\n({self.url})"
    
    def __eq__(self, value):
        if (self.name != value.name or
            self.url != value.url or 
            set(self.regions) != set(value.regions) or
            self.acronym != value.acronym):
            return False
        else:
            return True

@dataclass
class TaxidRegulation:
    """
    Container for regulatory information of a given taxid from
    a specific list. Structure data that is an expected output from
    the regulations module.
    """
    category : str = ""
    name : str = ""
    notes : str = ""
    derived_from : str = ""
    preferred_taxonomy_name : str = ""
    other_taxonomy_name : str = ""
    tax_id : int = 0
    list_acronym : str = ""
    target : str = ""
    hazard_group : str = ""
    uniprot : str = ""
    genbank_protein : str = ""

    @classmethod
    def from_row(cls, row: pd.Series, index_val=None):
        """
        Construct a TaxidRegulation from a pandas row returned by .iterrows().
        The optional index_val can be passed in if the DataFrame index should
        map to the `taxid` field.
        """
        # dataclass field names
        dataclass_fields = {f.name for f in fields(cls)}

        # build kwargs by selecting overlapping keys
        kwargs = {col: row[col] for col in row.index if col in dataclass_fields}

        # if the dataframe index is intended to represent `taxid`, override here
        if index_val is not None and "tax_id" in dataclass_fields:
            kwargs["tax_id"] = index_val

        return cls(**kwargs)

class RegulationType(StrEnum):
    """
    At what level a regulation is for i.e. "regulated at the species level."
    ??? These delineations require further thought before being practical.
    Specifically, the idea is to capture with an enum the type of regulation
    we are dealing with, which may derive some information about what it might have:
    taxid, uniprot, genbank accession. Whether or not it is pathogenic or if it targets
    something like Brazilian mushrooms.
    """
    INVASIVE = "invasive" # Non-pathogenic risk factors
    ORGANISM = "organism" # All species of this organism are regulated.
    SPECIES = "species" # This specific species is regulated.
    PEPTIDE = "peptide" # This toxin/protein/peptide is regulated, regardless of organism.

@dataclass
class RegulationOutput:
    """
    Container for regulation list information. Formatted for use in output JSON.

    * `name` ( str ) : Common name identification for this regulation list.
    * `region` ( str ) : 2 Letter Country codes for which this regulation applies
    * `regulation_level` ( str ) : At what level does this regulation apply to the hit (e.g. organism, species, protein)

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
REGULATED_TAXID_ANNOTATIONS : pd.DataFrame = pd.DataFrame({
    "category": pd.Series(dtype="str"),
    "name": pd.Series(dtype="str"),
    "notes": pd.Series(dtype="str"),
    "derived_from": pd.Series(dtype="str"),
    "preferred_taxonomy_name": pd.Series(dtype="str"),
    "other_taxonomy_name": pd.Series(dtype="str"),
    "tax_id": pd.Series(dtype="Int64"),
    "genbank_protein": pd.Series(dtype="str"),
    "uniprot": pd.Series(dtype="str"),
    "list_acronym": pd.Series(dtype="str"),
    "target": pd.Series(dtype="str"),
    "hazard_group": pd.Series(dtype="str")
    })

# Map of children taxids to the regulated taxid in the lists.
CHILD_TAXID_MAP = pd.DataFrame({
    "child_taxid": pd.Series(dtype="Int64"),
    "regulated_taxid": pd.Series(dtype="Int64")
    })

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
    expected_cols = set(REGULATED_TAXID_ANNOTATIONS.columns)
    essential_cols = {"list_acronym", "tax_id"}
    input_cols = set(input_data.columns)
    # Check for missing columns
    if not essential_cols.issubset(input_data.columns):
    #if not (essential_cols in input_cols):
        raise ValueError(f"Input data must contain columns {essential_cols}, "
                         f"got {set(input_data.columns)}")

    input_data = input_data.reindex(columns=expected_cols)
    # Coerce invalid entries ("TBD", empty, etc.) to NaN, then cast to Int64
    input_data["tax_id"] = (
        pd.to_numeric(input_data["tax_id"], errors="coerce")
        .astype("Int64")
    )

    # The input data requires dtypes to be specified for any missing info.
    input_data = input_data.astype(REGULATED_TAXID_ANNOTATIONS.dtypes.to_dict())
    REGULATED_TAXID_ANNOTATIONS = pd.concat(
        [REGULATED_TAXID_ANNOTATIONS, input_data],
        ignore_index=True)

def add_child_lut_data(input_data: pd.DataFrame):
    """
    Append validated child TaxID mappings to CHILD_TAXID_MAP.
    Ensures correct columns, restricts to them, and removes duplicates.
    """
    global CHILD_TAXID_MAP
    expected_cols = set(CHILD_TAXID_MAP.columns)

    # Check for missing columns
    if not expected_cols.issubset(input_data.columns):
        raise ValueError(f"Input data must contain columns {expected_cols}, "
                         f"got {list(input_data.columns)}")

    # Restrict to only the expected columns
    input_data = input_data.reindex(columns=expected_cols)
    input_data["child_taxid"] = pd.to_numeric(input_data["child_taxid"], errors="coerce").astype("Int64")
    input_data["regulated_taxid"] = pd.to_numeric(input_data["regulated_taxid"], errors="coerce").astype("Int64")
    input_data = input_data.astype(CHILD_TAXID_MAP.dtypes.to_dict())
    CHILD_TAXID_MAP = pd.concat(
        [CHILD_TAXID_MAP, input_data], ignore_index=True
        ).drop_duplicates().reset_index(drop=True)

def add_regulated_list(new_list : RegulationList) -> bool:
    """
    Wrapper for safely adding a list to global data.
    """
    global REGULATION_LISTS

    # Overwrite protection:
    existing = REGULATION_LISTS.get(new_list.acronym)
    # Save the list.
    if existing:
        if existing == new_list:
            REGULATION_LISTS[new_list.acronym] = new_list
            return True
        logger.warning("Attempting to assign different list information to the same"
        "list identifier:\nExisting:\n%sNew:\n%s. List will be skipped.", existing, new_list)
        return False

    REGULATION_LISTS[new_list.acronym] = new_list
    return True
