#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Set of containers for storing information important to screen
 outputs. Information is stored as a structure of dataclasses, and are
 converted between the dataclass / dict / json_file as required (using json.py module).

 The "annotations" dictionary, present in the HitResult,
 contains non-structured information, and is populated with differing information
 under differing keys depending on which step the information
 is derived (Biorisk, Taxonomy etc)

 In this way, the Results object serves as a common state, that can be updated
 whilst not being temporally appended like a log file i.e. .screen file.

 The Result all pertinent output information of a run.

 A Screen is made up of several Queries,
 which are made up of hits.
 Each hit derives from a Step, and has an associated recommendation.
 The Query's recommendation is the result of parsing all hitResults recommendations.

 ScreenResult:
     [QueryResult]
         recommendation (per query)
         [HitResult]:
             recommendation (per hit)
"""
import re
import logging
from dataclasses import dataclass, asdict, field
from typing import List, Iterator, Tuple
from enum import StrEnum
from importlib.metadata import version, PackageNotFoundError
import pandas as pd
from commec.tools.search_handler import SearchToolVersion

logger = logging.getLogger(__name__)

try:
    COMMEC_VERSION = version("commec")
except PackageNotFoundError:
    COMMEC_VERSION = "error"

# Seperate versioning for the output JSON.
JSON_COMMEC_FORMAT_VERSION = "0.1"


class ScreenStatus(StrEnum):
    """
    All possible outputs from commec for a query screen, or individual hit.
    Ordered by importance of user feedback, and by severity of screened outcome.
    """

    NULL = "-"
    SKIP = "Skip"
    PASS = "Pass"
    CLEARED_WARN = "Warning (Cleared)"
    CLEARED_FLAG = "Flag (Cleared)"
    WARN = "Warning"
    FLAG = "Flag"
    ERROR = "Error"

    @property
    def description(self) -> str:
        """Return a plaintext description of the status"""
        descriptions = {
            ScreenStatus.NULL: (
                "Status not initialized due to an error, interrupt,"
                " or other unexpected outcome"
            ),
            ScreenStatus.SKIP: (
                "Screening step intentionally skipped (e.g. skipping taxonomy screen in FAST mode,"
                " skipping benign screen when there are no flags to clear)"
            ),
            ScreenStatus.PASS: "Query was not flagged in this screening step; biosecurity review may not be needed",
            ScreenStatus.CLEARED_WARN: (
                "Warning was cleared, since query region was identified as benign"
                " (e.g. housekeeping gene, common synbio part)"),
            ScreenStatus.CLEARED_FLAG: (
                "Flag was cleared, since query region was identified as benign"
                " (e.g. housekeeping gene, common synbio part)"),
            ScreenStatus.WARN: (
                "Possible sequence of concern identified, but with low confidence"
                "(e.g. virulence factors or proteins shared among regulated and non-regulated organisms)"
            ),
            ScreenStatus.FLAG: "Query contains sequence of concern and requires additional biosecurity review",
            ScreenStatus.ERROR: "An error occured and this step failed to run",
        }
        return descriptions[self]

    @property
    def importance(self):
        """Encode the importance of each status."""
        order = {
            ScreenStatus.NULL: 0,
            ScreenStatus.SKIP: 1,
            ScreenStatus.PASS: 2,
            ScreenStatus.CLEARED_WARN: 3,
            ScreenStatus.CLEARED_FLAG: 4,
            ScreenStatus.WARN: 5,
            ScreenStatus.FLAG: 6,
            ScreenStatus.ERROR: 7,
        }
        return order[self]

    def clear(self):
        """Convert a WARN or FLAG into its cleared counterpart, and return that."""
        if self == ScreenStatus.WARN:
            return ScreenStatus.CLEARED_WARN
        if self == ScreenStatus.FLAG:
            return ScreenStatus.CLEARED_FLAG
        return self
    
    def __gt__(self, value):
        return self.importance > value.importance

    def __lt__(self, value):
        return self.importance < value.importance

    def __ge__(self, value):
        return self.importance >= value.importance

    def __le__(self, value):
        return self.importance <= value.importance

def compare(a: ScreenStatus, b: ScreenStatus):
    """
    Compare two recommendations, return the most important one.
    """
    if a.importance > b.importance:
        return a
    return b


class ScreenStep(StrEnum):
    """
    Enumeration of the Steps for Commec screening
    """

    BIORISK = "Biorisk Search"
    TAXONOMY_NT = "Nucleotide Taxonomy Search"
    TAXONOMY_AA = "Protein Taxonomy Search"
    BENIGN_PROTEIN = "Benign Protein Search"
    BENIGN_RNA = "Benign RNA Search"
    BENIGN_DNA = "Benign DNA Search"


@dataclass
class HitScreenStatus:
    """
    Maps ScreenStatus to ScreenStep for a single hit.
    """

    status: ScreenStatus = ScreenStatus.NULL
    from_step: ScreenStep = field(default_factory=ScreenStep)


@dataclass
class MatchRange:
    """
    Container for coordinate information of where hits match to a query.
    """

    e_value: float = 0.0
    # percent identity?
    match_start: int = 0
    match_end: int = 0
    query_start: int = 0
    query_end: int = 0

    # TODO: Add frame, as QueryStart and QueryEnd should be in Frame0 NT coords.

    def length(self):
        """
        Returns the length in Nucleotides of 
        this range for the query coordinates.
        """
        return abs(self.query_end - self.query_start)

    def __hash__(self):
        return hash(
            (
                self.e_value,
                self.match_start,
                self.match_end,
                self.query_start,
                self.query_end,
            )
        )

    def __eq__(self, other):
        if not isinstance(other, MatchRange):
            return NotImplemented
        return (
            self.e_value == other.e_value
            and self.match_start == other.match_start
            and self.match_end == other.match_end
            and self.query_start == other.query_start
            and self.query_end == other.query_end
        )
    
    def __str__(self):
        return f"{self.query_start}-{self.query_end}"


@dataclass
class HitResult:
    """
    Container for all information regarding a single hit with a range(s), to a single query.
    A Hit is any outcome from any step during Commec Screening, which can be mapped onto the Query.
    """

    recommendation: HitScreenStatus = field(default_factory=HitScreenStatus)
    name: str = ""
    description: str = ""
    ranges: list[MatchRange] = field(default_factory=list)
    annotations: dict = field(default_factory=dict)

    def get_e_value(self) -> float:
        """Gets the best e-value across all ranges, useful for sorting hits"""
        out: float = 10.0
        for r in self.ranges:
            out = min(out, r.e_value)
        return out

    def __str__(self) -> str:
        output = (
            f"{self.name}: {self.description}.\n{self.recommendation.from_step}"
            f", {self.recommendation.status}. Ranges (#{len(self.ranges)})\n"
            )
        match_string = ""
        for r in self.ranges:
            match_string += f"{r.query_start}-{r.query_end}, "
        match_string = match_string[:-2]

        return output + "[" + match_string + "]"

@dataclass
class QueryScreenStatus:
    """
    Summarises the status across all hits for a single query.
    """

    screen_status: ScreenStatus = ScreenStatus.NULL
    biorisk_status: ScreenStatus = ScreenStatus.NULL
    protein_taxonomy_status: ScreenStatus = ScreenStatus.NULL
    nucleotide_taxonomy_status: ScreenStatus = ScreenStatus.NULL
    benign_status: ScreenStatus = ScreenStatus.NULL

    def update_query_flag(self):
        """
        Parses the query status across all screening steps, then updates overall flag.
        Designed to be called at any time after Step 1.
        """
        if self.benign_status is not ScreenStatus.NULL:
            self.screen_status = self.benign_status
            return

        if self.biorisk_status == ScreenStatus.FLAG:
            self.screen_status = ScreenStatus.FLAG
            return

        if (
            self.protein_taxonomy_status == ScreenStatus.FLAG
            or self.nucleotide_taxonomy_status == ScreenStatus.FLAG
        ):
            self.screen_status = ScreenStatus.FLAG
            return

        if self.biorisk_status == ScreenStatus.WARN:
            self.screen_status = ScreenStatus.WARN
            return

        if (
            self.protein_taxonomy_status == ScreenStatus.WARN
            or self.nucleotide_taxonomy_status == ScreenStatus.WARN
        ):
            self.screen_status = ScreenStatus.WARN
            return

        # At the moment we only get here when biorisk screen is just run.
        # this will set the global recommend to pass or error, after biorisk is done.
        # We will need earlier checks to maintain the error status of future steps however.
        if self.screen_status == ScreenStatus.NULL:
            self.screen_status = self.biorisk_status


@dataclass
class QueryResult:
    """
    Container to hold screening result data pertinant to a single Query
    """

    query: str = ""
    length: int = 0
    #sequence: str = ""
    recommendation: QueryScreenStatus = field(default_factory=QueryScreenStatus)
    hits: dict[str, HitResult] = field(default_factory=dict)

    def get_hit(self, match_name: str) -> HitResult:
        """Wrapper for get logic."""
        return self.hits.get(match_name)
    
    def check_hit_range(self, input_region : MatchRange):
        """
        Checks all existing hits for whether there is an similar query coordinate region.
        Returns the relevant hit, or None.
        """
        for hit in self.hits.values():
            for region in hit.ranges:
                if (
                    input_region.query_start == region.query_start
                    and input_region.query_end == region.query_end
                ):
                    return hit
        return None

    def add_new_hit_information(self, new_hit: HitResult) -> bool:
        """
        Adds a Hit Description to this query, but only adds if the hit is unique, or has a new range.
        Returns True if the hit was not unique, but added unique info to the hit.
        """
        existing_hit = self.hits.get(new_hit.name)
        hits_is_updated: bool = False

        # Nothing matches in Name, try a matched region.
        if not existing_hit:
            for region in new_hit.ranges:
                existing_hit = self.check_hit_range(region)
                if existing_hit:
                    logger.debug("Using existing hit from shared region: %s", existing_hit)
                    hits_is_updated = True # We want to append info if new hit is differently named.
                    break

        # Nothing matches in Name or region... new hit!
        if not existing_hit:
            self.hits[new_hit.name] = new_hit
            return False

        for new_region in new_hit.ranges:
            is_unique_region = True
            for existing_region in existing_hit.ranges:
                if (
                    new_region.query_start == existing_region.query_start
                    and new_region.query_end == existing_region.query_end
                ):
                    is_unique_region = False
            if is_unique_region:
                hits_is_updated = True
                existing_hit.ranges.append(new_region)

        return hits_is_updated

    def get_flagged_hits(self) -> List[HitResult]:
        """
        Calculates and returns the list of hits, for all Warnings or Flags.
        Typically used as the regions to check against for benign screens.
        """
        flagged_and_warnings_data = [
            flagged_hit
            for flagged_hit in self.hits.values()
            if flagged_hit.recommendation.status
            in {ScreenStatus.WARN, ScreenStatus.FLAG}
        ]
        return flagged_and_warnings_data

    def _update_step_flags(self):
        """
        Updates the steps within QueryScreenStatus to be congruent for every hit.
        Then updates the Query flag to consolidate all.
        """
        for hit in self.hits.values():
            match hit.recommendation.from_step:
                case ScreenStep.BIORISK:
                    self.recommendation.biorisk_status = compare(
                        self.recommendation.biorisk_status, hit.recommendation.status)
                case ScreenStep.TAXONOMY_AA:
                    self.recommendation.protein_taxonomy_status = compare(
                        self.recommendation.protein_taxonomy_status, hit.recommendation.status)
                case ScreenStep.TAXONOMY_NT:
                    self.recommendation.nucleotide_taxonomy_status = compare(
                        self.recommendation.nucleotide_taxonomy_status, hit.recommendation.status)
                case ScreenStep.BENIGN_PROTEIN:
                    self.recommendation.benign_status = compare(
                        self.recommendation.benign_status, hit.recommendation.status)
                case ScreenStep.BENIGN_RNA:
                    self.recommendation.benign_status = compare(
                        self.recommendation.benign_status, hit.recommendation.status)
                case ScreenStep.BENIGN_DNA:
                    self.recommendation.benign_status = compare(
                        self.recommendation.benign_status, hit.recommendation.status)

        self.recommendation.update_query_flag()

    def update(self):
        """
        Call this before exporting to file.
        Ensures 
        Sorts the hits based on E-values,
        Updates the commec recommendation based on all hits recommendations.
        """
        # A rare instance where we want our dictionary to be sorted
        # So we convert it to a list, and recreate a dictionary - which maintains the order.
        # Thankfully we only do this once per Query at the end.
        sorted_items_desc = sorted(
            self.hits.items(), key=lambda item: item[1].get_e_value(), reverse=True
        )
        self.hits = dict(sorted_items_desc)
        self._update_step_flags()


@dataclass 
class SearchToolInfo:
    """ Container to hold version info for search tools and databases used. """
    biorisk_search_info:        SearchToolVersion = field(default_factory=SearchToolVersion)
    protein_search_info:        SearchToolVersion = field(default_factory=SearchToolVersion)
    nucleotide_search_info:     SearchToolVersion = field(default_factory=SearchToolVersion)
    benign_protein_search_info: SearchToolVersion = field(default_factory=SearchToolVersion)
    benign_rna_search_info:     SearchToolVersion = field(default_factory=SearchToolVersion)
    benign_dna_search_info:     SearchToolVersion = field(default_factory=SearchToolVersion)


@dataclass
class ScreenRunInfo:
    """Container dataclass to hold general run information for a commec screen"""
    commec_version: str = str(COMMEC_VERSION)
    json_output_version: str = JSON_COMMEC_FORMAT_VERSION
    time_taken: str = ""
    date_run: str = ""
    search_tool_info: SearchToolInfo = field(default_factory=SearchToolInfo)

@dataclass
class ScreenQueryInfo:
    """ Container for summarising the query input data """
    file: str = ""
    number_of_queries: int = 0
    total_query_length: int = 0

@dataclass
class ScreenResult:
    """
    Root dataclass to hold all data related to the screening of an individual query by commec.
    """
    commec_info: ScreenRunInfo = field(default_factory=ScreenRunInfo)
    query_info: ScreenQueryInfo = field(default_factory=ScreenQueryInfo)
    queries: dict[str, QueryResult] = field(default_factory=dict)

    def get_query(self, query_name: str) -> QueryResult:
        """
        Wrapper for Query get logic.
        """
        search_term = query_name
        if re.search(r'_[1-6]$', query_name):  # Check if string ends with _1 to _6
            search_term = query_name[:-2]  # Remove last two characters

        return self.queries.get(search_term)

    def update(self):
        """
        Propagate update to all children dataclasses.
        """
        for query in self.queries.values():
            query.update()

    def regions(self) -> Iterator[Tuple[QueryResult, HitResult, MatchRange]]:
        """
        Helper function, iterates through all queries, hits, and regions in the ScreenResult object.
        Yields tuples of (query, hit, region).
        """
        for query in self.queries.values():
            for hit in query.hits.values():
                for region in hit.ranges:
                    yield query, hit, region

    def hits(self) -> Iterator[Tuple[QueryResult, HitResult]]:
        """
        Helper function, iterates through all queries and hits in the ScreenResult object.
        Yields tuples of (query, hit).
        """
        for query in self.queries.values():
            for hit in query.hits.values():
                yield query, hit

    def get_flag_data(self) -> pd.DataFrame:
        """
        Returns a dataframe containing the status' from each screen step.
        Useful for printing summary information.
        """
        data = []
        for query in self.queries.values():
            data.append({
                "query": query.query[:25],
                "overall": query.recommendation.screen_status,
                "biorisk": query.recommendation.biorisk_status,
                "taxonomy_aa": query.recommendation.protein_taxonomy_status,
                "taxonomy_nt": query.recommendation.nucleotide_taxonomy_status,
                "benign": query.recommendation.benign_status
            })
        
        output_data : pd.DataFrame = pd.DataFrame(data)
        return output_data

    def __str__(self):
        df = self.get_flag_data()
        return df.to_string(index=False, col_space = 12)
    
    def __repr__(self):
        return str(asdict(self))
