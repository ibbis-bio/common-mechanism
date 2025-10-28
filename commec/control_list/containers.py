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
import re
from dataclasses import dataclass, field, fields
from typing import Optional
from enum import StrEnum
import pandas as pd

logger = logging.getLogger(__name__)

class AccessionFormat(StrEnum):
    """
    Supported hashable accesion formats.
    """
    TAXID = "TaxonomyID"
    GENBANK = "GenbankProtein"
    UNIPROT = "UniprotID"

    def __repr__(self):
        return f"<Fmt:{self}>"

def derive_accession_type(identifier : str) -> AccessionFormat | None:
    """
    Utilises regex to determine whether the input
    string is a TaxID, Uniprot, or Genbank accesion.
    Returns None if it cannot be determined.
    """
    if isinstance(identifier, int):
        return AccessionFormat.TAXID

    patterns = {
        AccessionFormat.UNIPROT: re.compile(
            r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
        ),
        AccessionFormat.TAXID: re.compile(
            r"^[0-9]+$"
        ),
        AccessionFormat.GENBANK: re.compile(
            r"^(?:[A-Z][0-9]{5}|[A-Z]{2}[0-9]{6,8}|[A-Z]{3}[0-9]{5,7})(?:\.\d+)?$"
        ),
    }

    for fmt, pattern in patterns.items():
        if pattern.fullmatch(identifier):
            return fmt

    return None

@dataclass(frozen=True)
class Accession:
    """
    A unified hashable identifier for 
    taxonomy IDs, GenBank, or UniProt accessions.

    Used for hashed index accession for rapid annotations retrieval
    from a pandas DataFrame.
    """
    code: str
    type: AccessionFormat   # one of: "taxid", "genbank", "uniprot"

    def __init__(self, accession: Optional[str]=None,
                       taxid: Optional[int]=None,
                       uniprot: Optional[str]=None,
                       genbank: Optional[str]=None,
                       row : pd.Series = None):
        # Best guess if unknown input.
        if accession:
            acc_fmt = derive_accession_type(accession)
            if acc_fmt:
                object.__setattr__(self, "code", accession)
                object.__setattr__(self, "type", acc_fmt)
                return

        # ensure only one is set
        values = [(taxid, AccessionFormat.TAXID),
                  (uniprot, AccessionFormat.UNIPROT),
                  (genbank, AccessionFormat.GENBANK)]
        non_empty = [(v, t) for v, t in values if not pd.isna(v)]
        if len(non_empty) == 0:
            error_id = row.to_string() if not row is None else "[Unknown row info]"
            logger.debug("Accession [\n%s\n] must have at least one non-empty"
                         " value for tax_id, genbank, or uniprot.\n"
                         "It has been assigned to TaxID 0.", error_id)
            non_empty = [("0", AccessionFormat.TAXID)]

        object.__setattr__(self, "code", str(non_empty[0][0]))
        object.__setattr__(self, "type", non_empty[0][1])

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
