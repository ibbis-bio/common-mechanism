# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

"""
Defines Containers to show default layouts for regulation data.
What is the shape of the regulation data, and where to store it.

    RegulationContainer : 
        A standard output that may be associated with a hit.
        Contains list, region, and level information.
    RegulationLevel :
        Determines valid level information strings.
    ControlListInfoContainer :
        Storage of all regulation data for a given taxid from a given list.

Also defines the storage data for 
    * regulation lists, 
    * per taxid regulation information per list, 
    * and child to parent taxid mapping.
"""
import re
from dataclasses import dataclass, field, fields
from typing import Optional
from enum import StrEnum
import pandas as pd

class AccessionFormat(StrEnum):
    """
    Supported hashable accesion formats.
    """
    TAXID = "TaxonomyID"
    NOT_SET = "-"
    UNKNOWN = "Unknown"

    def __repr__(self):
        return f"<Fmt:{self}>"


TAXID_PATTERN = re.compile(r"^[0-9]+$")

def derive_accession_format(accession: str) -> AccessionFormat:
    """
    Determine supported accession formats using regex-based matching.

    Returns:
        AccessionFormat.TAXID if numeric (ASCII digits only),
        AccessionFormat.UNKNOWN otherwise.
    """
    if not accession:
        return AccessionFormat.UNKNOWN

    accession_str = str(accession)

    patterns = {
        AccessionFormat.TAXID: TAXID_PATTERN,
    }

    # Find and return the first matching format, or UNKNOWN if none match
    return next(
        (fmt for fmt, pattern in patterns.items() if pattern.fullmatch(accession_str)),
        AccessionFormat.UNKNOWN
    )


@dataclass(frozen=True)
class Accession:
    """
    A unified hashable wrapper for arbitrary accession types. e.g. TaxID

    Used for hashed index accession for rapid annotations retrieval
    from a pandas DataFrame. 
    
    Includes helper functions to pre-calculate
    what type of accession this is for user reporting.
    """
    code : str
    type : AccessionFormat

    def __init__(self, accession: Optional[str]=None):
        """
        An Accession is created from an arbitrary str-like with truthiness.
        Invalid accessions will instantiated as None.
        """
        if accession:
            object.__setattr__(self, "code", str(accession))
            object.__setattr__(self, "type", derive_accession_format(accession))
            return

        object.__setattr__(self, "code", None)
        object.__setattr__(self, "type", AccessionFormat.UNKNOWN)

    def get_format(self) -> AccessionFormat:
        """
        Returns the likely format of this accession,
        calculates it if it hasn't been calculated yet.
        """
        return self.type

    def __hash__(self):
        return hash(self.code)

    def get(self) -> tuple[str, str]:
        """Return (code, type)."""
        return self.code, self.type

    def __str__(self):
        return f"{self.code}"

    def __repr__(self):
        return f"{self.type}:{self.code}"

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
    
    def __hash__(self):
        return hash(self.acronym)

class ListMode(StrEnum):
    """
    Describes how a list is interpreted.
    When an Accession is identified as part of a list, whether or not we treat
    that taxid as regulated or not will depend on the following modes:
    
    * COMPLIANCE - All taxids from this list are regulated.
    * CONDITIONAL_NUM - Only mark taxid as regulated if it also appears in another list.
    * COMPLIANCE_WARN - Mark Taxid as only a WARNING.
    * IGNORE - Ignore this list entirely.

    Future settings can be used to apply these to the lists under various
    conditions, however for now, we simply set any list that is not
    part of the regional context, to CONDITIONAL_NUM, as default Commec
    behaviour.
    """
    COMPLIANCE = "Compliance"
    CONDITIONAL_NUM = "Conditional Compliance"
    COMPLIANCE_WARN = "Comply with Warning"
    IGNORE = "Ignored"

@dataclass
class ControlList:
    """
    Contains the name, acronym, url, and affected regions for a Control List.
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
class ControlListInfo:
    """
    Container for regulatory information of a given accession from
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
        Construct a ControlListInfo from a pandas row returned by .iterrows().
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

class CategoryType(StrEnum):
    """
    At what level a regulation is for i.e. "regulated at the species level."
    ??? These delineations require further thought before being practical.
    Specifically, the idea is to capture with an enum the type of regulation
    we are dealing with, which may derive some information about what it might have:
    taxid, uniprot, genbank accession. Whether or not it is pathogenic or if it targets
    something like Brazilian mushrooms.
    This will be used to parse the category column of the input csvs.
    """
    BACTERIA = "Bacteria"
    VIRUSES = "Viruses"
    EUKARYOTA = "Eukaryota"
    PROTEIN = "Proteins"
    TOXIN = "Toxins"
    OTHER_EUKARYOTA = "Other Eukaryota"
    OTHER_EUKAYROTA_ANIMAL = "Other Eukaryota (Animal)"
    NON_PROTEIN_TOXIN = "Non-Protein Toxin"
    TOXIN_SYNTHESIS_ENZYME = "Toxin Synthesis Enzyme"
    PRIONS_AND_TSE = "Prions & TSEs"
    NONE = "None"

@dataclass
class RegulationOutput:
    """
    Container for regulation list information. Formatted for use in output JSON.

    * `list` ( str ) : Common acroynm for control list.
    * `category` ( str ) : What category does this control list apply to the hit (e.g. organism, species, protein)

    See results.py for other examples.
    """
    name : str = ""
    category : str = ""
    # regulation_level : RegulationLevel = field(default_factory=RegulationLevel) # At what level this list is regulated.
    # toxicity : str = "no data" - TBD
