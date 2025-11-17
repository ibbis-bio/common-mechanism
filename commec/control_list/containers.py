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

class ListUseAcronym(StrEnum):
    """
    Storage for the acroynms of the various formats of ControlList types.
    """
    EXPORTCONTROLS	= "EXPORT"
    FRAMEWORK	  	= "FRMWK"
    DUALUSE	    	= "DURC"
    GMO			    = "GMO"
    OTHERPATHOGEN	= "PATHGN"
    BENIGN  		= "BENIGN"

@dataclass
class ControlList:
    """
    Contains the name, acronym, url, and affected regions for a Control List.
    """
    name : str = ""
    acronym : str = ""
    url : str = ""
    region : Region = field(default_factory=Region)
    status : ListMode = field(default_factory=ListMode)
    use : ListUseAcronym = field(default_factory=ListUseAcronym)

    def __str__(self):
        """
        Shorthand acroynms for list comprehension. e.g:
        AG_CCL_EXPORT, or IN_SCOMET_EXPORT.
        Describes the region affect, the acronym, and the use in a single line.
        For JSON and logging reporting purposes.
        """
        return f"{self.region.acronym}_{self.acronym}_{self.use}"
    
    def description(self) -> str:
        """
        Human readable text based description of this control list.
        """
        return f"[{self.acronym}] {self.name} - {self.region.name}\n({self.url})"
    
    def __eq__(self, value):
        if (self.name != value.name or
            self.url != value.url or
            self.region != value.region or
            self.acronym != value.acronym):
            return False
        else:
            return True

class CategoryType(StrEnum):
    """
    Valid options for values under the 'category' column in regulated_taxids.csv
    inputs. Communicates the type of entity being references by the control list.
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
    NONE = "None" # Default value.

@dataclass
class ControlListContext:
    """
    Provides additional context for how the results of querying a control list
    applies to the queried accession.

    *`derived_from` (str) : The name of the entity referenced by the Control List, if different from the taxonomy name.
    *`is_child` (bool) : Whether or not the queried taxid was a child of the controlled taxid.
    """
    derived_from : str = None
    is_child : bool = False
    #should_ignore : bool = False

@dataclass
class ControlListOutput:
    """
    Container for Control List information. Formatted for use in output JSON.

    * `list` ( str ) : Common acroynm for control list.
    * `name` ( str ) : High level name used for this annotation.
    * `category` ( str ) : What category does this control list apply to the hit (e.g. organism, species, protein)

    See results.py for other examples.
    """
    name : str = ""
    category : str = ""
    list : str = ""
