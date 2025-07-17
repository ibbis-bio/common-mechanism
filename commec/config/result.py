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
from commec import __version__ as COMMEC_VERSION

logger = logging.getLogger(__name__)



# Seperate versioning for the output JSON.
JSON_COMMEC_FORMAT_VERSION = "0.3"


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
                " skipping low-concern screen when there are no flags to clear)"
            ),
            ScreenStatus.PASS: "Query was not flagged in this screening step; biosecurity review may not be needed",
            ScreenStatus.CLEARED_WARN: (
                "Warning was cleared, since query region was identified as low-concern"
                " (e.g. housekeeping gene, common synbio part)"),
            ScreenStatus.CLEARED_FLAG: (
                "Flag was cleared, since query region was identified as low-concern"
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
    LOW_CONCERN_PROTEIN = "Benign Protein Search"
    LOW_CONCERN_RNA = "Benign RNA Search"
    LOW_CONCERN_DNA = "Benign DNA Search"


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

class Rationale(StrEnum):
    """
    Container for rationale texts in Commec outputs in one place.

    When reporting rationale, the primary (see ScreenStatus.importance attribute)
    status is always reported first. After this, secondary less important statuses are
    reported (usually warnings) in the rationale.
    """
    NULL = "-"
    ERROR = "There was an error during "

    # Start rationale
    START_PRIMARY = "Matches "
    START_PASS = "No regions of concern"
    START_SECONDARY = "; as well as "

    # pre types:
    BIORISK_FLAG = "pathogenic or toxin function"
    BIORISK_WARN = "virulence factor"

    # Taxonomy types
    PR = "protein"
    NT = "nucleotide"

    BODY = "sequence with"

    # post Types:
    TAX_FLAG = " regulated organisms"
    TAX_WARN = " equally-good matches to regulated and non-regulated organisms"


    # Outcomes:
    NOTHING = ("No matches found during any stage of analysis. "
                "Sequence risk is unknown, possibly generated in silico. ")
    TOO_SHORT = "Query is too short, and was skipped."

    FLAG = " flags"
    WARN = " warnings"
    FLAGWARN = " flags and warnings"
    CLEARED = " cleared as common or non-hazardous"


@dataclass
class QueryScreenStatus:
    """
    Summarises the status across all hits for a single query.
    """

    screen_status: ScreenStatus = ScreenStatus.NULL
    biorisk: ScreenStatus = ScreenStatus.NULL
    protein_taxonomy: ScreenStatus = ScreenStatus.NULL
    nucleotide_taxonomy: ScreenStatus = ScreenStatus.NULL
    low_concern: ScreenStatus = ScreenStatus.NULL
    rationale : str = "-"

    # Mapping between screen steps and the fields above
    STEP_TO_STATUS_FIELD = {
        ScreenStep.BIORISK: 'biorisk',
        ScreenStep.TAXONOMY_NT: 'nucleotide_taxonomy', 
        ScreenStep.TAXONOMY_AA: 'protein_taxonomy',
        ScreenStep.LOW_CONCERN_PROTEIN: 'low_concern',
        ScreenStep.LOW_CONCERN_RNA: 'low_concern',
        ScreenStep.LOW_CONCERN_DNA: 'low_concern',
    }

    def update_step_status(self, step: ScreenStep, status: ScreenStatus, override_skip: bool = False) -> None:
        """
        Update the query screen status for a particular step if the proposed status is more
        important than the current one.
        """
        _, current_status = self._get_step_field_and_status(step)
        if status.importance > current_status.importance:
            self.set_step_status(step, status, override_skip)

    def set_step_status(self, step: ScreenStep, status: ScreenStatus, override_skip: bool = False) -> None:
        """
        Set the query screen status for a particular step.
        In most cases, query steps that have already been skipped should not be updated.
        """
        field_name, current_status = self._get_step_field_and_status(step)
        if override_skip or current_status != ScreenStatus.SKIP:
            setattr(self, field_name, status)

    def _get_step_field_and_status(self, step: ScreenStep) -> tuple[str, ScreenStatus]:
        field_name = QueryScreenStatus.STEP_TO_STATUS_FIELD.get(step)
        return field_name, getattr(self, field_name)

    def update(self, query_data):
        """
        Updates the overall status flag for this query, based on the status
        from each step. Some special cases are also handled:
        * Skipping this query entirely
        * This query passed, but is suspiciously new.
        ----
        Inputs:
        query_data : Query - The input Query as loaded by Screen, see Query.py
        """
        self.screen_status = max(
            self.biorisk,
            self.protein_taxonomy,
            self.nucleotide_taxonomy,
            self.low_concern
        )

        # If everything is happy, but we haven't hit anything, time to be suspicious...
        if (self.screen_status == ScreenStatus.PASS and query_data.no_hits_warning):
            self.screen_status = ScreenStatus.WARN
            return

        # If biorisk was skipped then it is skipped overall - likely query is too short...
        if (self.biorisk == ScreenStatus.SKIP):
            self.screen_status = ScreenStatus.SKIP
            return


    def __str__(self) -> str:
        output = f"""
                Overall     : {self.screen_status}\n
                Biorisk     : {self.biorisk}\n
                Protein     : {self.protein_taxonomy}\n
                Nucleotide  : {self.nucleotide_taxonomy}\n
                Low Concern : {self.low_concern}\n
                {self.rationale}
                """
        return output

    def get_error_stepname(self):
        """
        Returns a text step name of the first error occurance for use in logging.
        """
        if self.biorisk == ScreenStatus.ERROR:
            return "Biorisk Screening"
        if self.protein_taxonomy == ScreenStatus.ERROR:
            return "Protein Taxonomy Screening"
        if self.nucleotide_taxonomy == ScreenStatus.ERROR:
            return "Nucleotide Taxonomy Screening"
        if self.low_concern == ScreenStatus.ERROR:
            return "Low concern Screening"

        return "Screening" # General Error at some stage.


@dataclass
class QueryResult:
    """
    Container to hold screening result data pertinant to a single Query
    """
    query: str = ""
    length: int = 0
    status: QueryScreenStatus = field(default_factory=QueryScreenStatus)
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
        Typically used as the regions to check against for low-concern screens.
        """
        flagged_and_warnings_data = [
            flagged_hit
            for flagged_hit in self.hits.values()
            if flagged_hit.recommendation.status
            in {ScreenStatus.WARN, ScreenStatus.FLAG}
        ]
        return flagged_and_warnings_data

    def _update_step_flags(self, query_data):
        """
        Updates the steps within QueryScreenStatus to be congruent for every hit.
        Then updates the Query flag to consolidate all.
        """
        logger.debug("Updating step status flags for query %s", self.query)
        logger.debug("Current status %s", self.status)
        
        ignored_status = {ScreenStatus.SKIP, ScreenStatus.ERROR, ScreenStatus.PASS}
        
        if self.status.biorisk not in ignored_status:
            self.status.biorisk = ScreenStatus.NULL
        if self.status.protein_taxonomy not in ignored_status:
            self.status.protein_taxonomy = ScreenStatus.NULL
        if self.status.nucleotide_taxonomy not in ignored_status:
            self.status.nucleotide_taxonomy = ScreenStatus.NULL
        if self.status.low_concern not in ignored_status:
            self.status.low_concern = ScreenStatus.NULL

        # Track status sets for rationale (only for specific steps that need it)
        status_sets = {
            ScreenStep.BIORISK: set(),
            ScreenStep.TAXONOMY_AA: set(), 
            ScreenStep.TAXONOMY_NT: set(),
        }

        # Collapse data from all hits:
        for hit in self.hits.values():
            step = hit.recommendation.from_step
            hit_status = hit.recommendation.status
        
            self.status.update_step_status(step, hit_status, override_skip=True)

            if step in status_sets:
                status_sets[step].add(hit_status)

        # Update Benign outcome based on the worst step.
        self.status.low_concern = max(
            self.status.low_concern,
            self.status.biorisk,
            self.status.protein_taxonomy,
            self.status.nucleotide_taxonomy
        )

        self.status.update(query_data)
        self._update_rationale(status_sets[ScreenStep.BIORISK],
                               status_sets[ScreenStep.TAXONOMY_AA],
                               status_sets[ScreenStep.TAXONOMY_NT])

        logger.debug("Updated status %s", self.status)

    def _update_rationale(self,
                          biorisks : set[ScreenStatus],
                          tax_aa : set[ScreenStatus],
                          tax_nt : set[ScreenStatus]):
        """ 
        Check existing statuses, and updates rationale accordingly.
        Requires sets containing unique statuses from each step, as
        each step is the primary status only. Passing all options 
        allows for more depth in rationale texts.
        """

        logger.debug("Biorisk set: %s", biorisks)
        logger.debug("TAX AA set: %s", tax_aa)
        logger.debug("TAX NT set: %s", tax_nt)

        state = self.status # Shorthand, accessor to be updated
        tax_all = tax_aa | tax_nt # Check both Taxonomy steps at once

        has_flags = state.screen_status == ScreenStatus.FLAG
        has_warns = ScreenStatus.WARN in biorisks | tax_aa | tax_nt
        has_clears = (ScreenStatus.CLEARED_FLAG in tax_all or
                      ScreenStatus.CLEARED_WARN in tax_all)

        logger.debug("%s has flags [%s], and has warnings [%s], and has clears [%s]",
                     self.query, has_flags, has_warns, has_clears)

        if state.screen_status == ScreenStatus.ERROR:
            state.rationale = Rationale.ERROR + state.get_error_stepname()
            return

        if state.screen_status == ScreenStatus.SKIP:
            state.rationale = Rationale.TOO_SHORT
            return

        # Unique circumstance, where there are no hits at all.
        if (state.screen_status == ScreenStatus.WARN and
            state.biorisk == ScreenStatus.PASS and
            state.protein_taxonomy == ScreenStatus.PASS and
            state.nucleotide_taxonomy == ScreenStatus.PASS and
            state.low_concern == ScreenStatus.PASS):
            state.rationale = Rationale.NOTHING
            return

        # Handle simple passes
        if state.screen_status == ScreenStatus.PASS:
            state.rationale = Rationale.START_PASS + "."
            return

        # Calculate any cleared outputs:
        tokens = ""
        if ScreenStatus.CLEARED_FLAG in tax_all:
            tokens = Rationale.FLAG
        if ScreenStatus.CLEARED_WARN in tax_all:
            tokens = Rationale.WARN
        if (ScreenStatus.CLEARED_FLAG in tax_all and
            ScreenStatus.CLEARED_WARN in tax_all):
            tokens = Rationale.FLAGWARN
        cleared_sentence = tokens + Rationale.CLEARED

        # Handle ONLY cleared outputs
        if state.screen_status in [ScreenStatus.CLEARED_FLAG,
                                   ScreenStatus.CLEARED_WARN]:
            state.rationale = Rationale.START_PASS + cleared_sentence
            return
        
        # Handle complex outputs:

        # Start creating rationale message:
        output = Rationale.START_PRIMARY
    
        types = []
        tax_types = []

        if has_flags:
            # "Matches FLAGS as well as WARNS"
            prebody = ""

            if ScreenStatus.FLAG in biorisks:
                logger.debug("Adding Biorisk Flag to primary.")
                types.append(Rationale.BIORISK_FLAG)
                prebody = Rationale.BODY + " "

            if ScreenStatus.FLAG in tax_aa:
                tax_types.append(Rationale.PR)
            if ScreenStatus.FLAG in tax_nt:
                tax_types.append(Rationale.NT)
            tax_types = " and ".join(tax_types)

            if ScreenStatus.FLAG in tax_all:
                types.append(tax_types + " " + Rationale.BODY + Rationale.TAX_FLAG)
            
            output += prebody + oxford_comma(types)

            if has_warns:
                output += Rationale.START_SECONDARY

        types = []
        tax_types = []

        if has_warns:
            # "Matches WARNS"
            prebody = ""
            if ScreenStatus.WARN in biorisks:
                types.append(Rationale.BIORISK_WARN)
                prebody = Rationale.BODY + " "


            if ScreenStatus.WARN in tax_aa:
                tax_types.append(Rationale.PR)
            if ScreenStatus.WARN in tax_nt:
                tax_types.append(Rationale.NT)
            tax_types = " and ".join(tax_types)

            if ScreenStatus.WARN in tax_all:
                types.append(tax_types + " " + Rationale.BODY + Rationale.TAX_WARN)

            if has_flags:
                prebody = ""

            output += prebody + oxford_comma(types)
        
        if has_clears:
            output += Rationale.START_SECONDARY[:-1] + cleared_sentence

        state.rationale = output + "."
        return

    
    def update(self, query_data):
        """
        Call this before exporting to file.
        Ensures 
        Sorts the hits based on E-values,
        Updates the commec recommendation based on all hits recommendations.
        """
        
        assert hasattr(query_data, "no_hits_warning")

        # A rare instance where we want our dictionary to be sorted
        sorted_items_desc = sorted(
            self.hits.items(), key=lambda item: item[1].get_e_value(), reverse=True
        )
        self.hits = dict(sorted_items_desc)
        self._update_step_flags(query_data)

    def skip(self):
        """
        Called to skip this query, sets all recommendations to skip.
        """
        self.status.screen_status = ScreenStatus.SKIP
        self.status.biorisk = ScreenStatus.SKIP
        self.status.protein_taxonomy = ScreenStatus.SKIP
        self.status.nucleotide_taxonomy = ScreenStatus.SKIP
        self.status.low_concern = ScreenStatus.SKIP
        logger.debug("Query %s has all statuses assigned to SKIP.", self.query)


@dataclass
class SearchToolInfo:
    """ Container to hold version info for search tools and databases used. """
    biorisk_search_info:        SearchToolVersion = field(default_factory=SearchToolVersion)
    protein_search_info:        SearchToolVersion = field(default_factory=SearchToolVersion)
    nucleotide_search_info:     SearchToolVersion = field(default_factory=SearchToolVersion)
    low_concern_protein_search_info: SearchToolVersion = field(default_factory=SearchToolVersion)
    low_concern_rna_search_info:     SearchToolVersion = field(default_factory=SearchToolVersion)
    low_concern_dna_search_info:     SearchToolVersion = field(default_factory=SearchToolVersion)


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

    def update(self, queries_data):
        """
        Propagate update to all children dataclasses.
        """
        for query_name, query in self.queries.items():
            query.update(queries_data[query_name])

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
                "overall": query.status.screen_status,
                "biorisk": query.status.biorisk,
                "taxonomy_aa": query.status.protein_taxonomy,
                "taxonomy_nt": query.status.nucleotide_taxonomy,
                "cleared": query.status.low_concern
            })

        output_data : pd.DataFrame = pd.DataFrame(data)
        return output_data

    def get_rationale_data(self) -> pd.DataFrame:
        """
        Returns a dataframe containing the overall statuse across steps,
        as well as a human readable rationale text.
        Useful for printing summary information.
        """
        data = []
        for query in self.queries.values():
            data.append({
                "query": query.query[:25],
                "overall": query.status.screen_status,
                "rationale": query.status.rationale,
            })

        output_data : pd.DataFrame = pd.DataFrame(data)
        return output_data
    
    def rationale_text(self) -> str:
        """ Outputs the rationale data as formatted text. """
        output = ""
        for row in self.get_rationale_data().itertuples(index=False):
            output += f"{row.query:<26}: {row.overall:<12} --> {row.rationale}\n"
        return output

    def flag_text(self) -> str:
        """ Outputs the flag table data as formatted text."""
        return self.get_flag_data().to_string(index=False, col_space = 12, line_width=2048)

    def __str__(self):
        return self.flag_text()

    def __repr__(self):
        return str(asdict(self))

def oxford_comma(inputs : list[str]) -> str:
    """
    Takes a list of strings: 
        * `[a,b,c]`, 
    and outputs a single formatted string:
        * `\"a, b, and c\"`
    """
    if len(inputs) == 0:
        return ""
    if len(inputs) == 1:
        return inputs[0]
    output = ""
    for text in inputs[:-1]:
        output += text + ", "
    return output + "and " + inputs[-1]
