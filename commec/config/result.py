#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
'''
   Set of containers for storing information important to screen
    outputs. Information is stored as a structure of dataclasses, and are
    converted between the dataclass / dict / json_file as required (using json.py module). 

    The "annotations" dictionary, present in the HitDescription, 
    contains non-structured information, and is populated with differing information 
    under differing keys depending on which step the information
    is derived (Biorisk, Taxonomy etc)

    In this way, the Results object serves as a common state, that can be updated
    whilst not being temporally appended like a log file i.e. .screen file.

    The Result all pertinent output information of a run.
'''

from dataclasses import dataclass, asdict, field
from typing import List, Iterator, Tuple
from enum import StrEnum
from commec.tools.search_handler import SearchToolVersion

# Seperate versioning for the output JSON.
JSON_COMMEC_FORMAT_VERSION = "0.1"

class CommecScreenStep(StrEnum):
    """
    Enumeration of possible stages of Commec screening
    """
    BIORISK = 'Biorisk Screen'
    TAXONOMY_NT = "Nucleotide Taxonomy Screen"
    TAXONOMY_AA = "Protein Taxonomy Screen"
    BENIGN_PROTEIN = "Benign Protein Screen"
    BENIGN_RNA = "Benign RNA Screen"
    BENIGN_SYNBIO = "Benign SynBio Screen"

class CommecRecommendation(StrEnum):
    """
    All possible outputs from commec for a query or individual hit.
    Ordered by importance of user feedback, and by severity of screened outcome.
    """
    NULL = '-' # This was not set.
    SKIP = 'Skip' # Intentionally skipped this step.
    PASS = 'Pass' # Commec has approved this query at this step.
    CLEARED_WARN = 'Warn (Cleared)' # Commec has cleared this warning during Benign Screen
    CLEARED_FLAG = 'Flag (Cleared)' # Commec has cleared this flag during Benign Screen
    WARN = 'Warn' # This query may be suspicious...
    FLAG = 'Flag' # Commec has flagged this query.
    ERROR = 'Error' # An error occured, such that this step failed to run.

    @property
    def importance(self):
        """ Encode the importance of each Recommendations value. """
        order = {
            CommecRecommendation.NULL: 0,
            CommecRecommendation.SKIP: 1,
            CommecRecommendation.PASS: 2,
            CommecRecommendation.CLEARED_WARN: 3,
            CommecRecommendation.CLEARED_FLAG: 4,
            CommecRecommendation.WARN: 5,
            CommecRecommendation.FLAG: 6,
            CommecRecommendation.ERROR: 7,
        }
        return order[self]

    def clear(self):
        """ Convert a WARN or FLAG into its cleared counterpart, and return that. """
        if self == CommecRecommendation.WARN:
            return CommecRecommendation.CLEARED_WARN
        if self == CommecRecommendation.FLAG:
            return CommecRecommendation.CLEARED_FLAG
        return self

def compare(a : CommecRecommendation, b : CommecRecommendation):
    """ Compare two recommendations, return the most important one. """
    if a.importance > b.importance:
        return a
    return b

@dataclass
class CommecScreenStepRecommendation:
    """ Pairs a recommendation with a Screening step."""
    outcome : CommecRecommendation = CommecRecommendation.NULL
    from_step : CommecScreenStep = field(default_factory=CommecScreenStep)

@dataclass
class CommecSummaryStatistics:
    """
    Useful data on the summary of a Commec Screen.
    """
    virus_hits : int = 0
    bacteria_hits : int = 0
    eukaryote_hits : int = 0

@dataclass
class MatchRange:
    """
    Container for information on matching a hit to a query
    """
    e_value : float = 0.0
    # percent identity?
    match_start : int = 0
    match_end : int = 0
    query_start : int = 0
    query_end : int = 0

    def query_length(self):
        return abs(self.query_start - self.query_end)

    def __hash__(self):
        return hash((self.e_value,
                     self.match_start,
                     self.match_end,
                     self.query_start,
                     self.query_end))

    def __eq__(self, other):
        if not isinstance(other, MatchRange):
            return NotImplemented
        return (self.e_value == other.e_value and
                self.match_start == other.match_start and
                self.match_end == other.match_end and
                self.query_start == other.query_start and
                self.query_end == other.query_end)

@dataclass
class HitDescription:
    """ Container for all information regarding a single hit with a range(s), to a single query."""
    recommendation : CommecScreenStepRecommendation = field(default_factory=CommecScreenStepRecommendation)
    name : str = ""
    description : str = ""
    ranges : list[MatchRange] = field(default_factory=list)
    annotations : dict = field(default_factory = dict)

    def get_e_value(self) -> float:
        """ Gets the best e-value across all ranges, useful for sorting hits"""
        out : float = 10.0
        for r in self.ranges:
            out = min(out, r.e_value)
        return out
    
@dataclass
class CommecRecommendationContainer:
    """
    Summarises the recommendations across all hits for a single query.
    """
    commec_recommendation : CommecRecommendation = CommecRecommendation.NULL
    biorisk_screen : CommecRecommendation = CommecRecommendation.NULL
    protein_taxonomy_screen : CommecRecommendation = CommecRecommendation.NULL
    nucleotide_taxonomy_screen : CommecRecommendation = CommecRecommendation.NULL
    benign_screen : CommecRecommendation = CommecRecommendation.NULL

    def update_commec_recommendation(self, list_of_hits : list[HitDescription]):
        """
        Parses the current state of the json,
        and updates the global commec recommendation.

        TODO: Consider the value of iterating over all HitDescriptions (provided as input)
        and updating things with stricter logic.
        """
        if self.benign_screen in {CommecRecommendation.CLEARED_FLAG,
                                  CommecRecommendation.CLEARED_WARN,
                                  CommecRecommendation.PASS}:
            self.commec_recommendation = self.benign_screen
            return

        if self.biorisk_screen == CommecRecommendation.FLAG:
            self.commec_recommendation = CommecRecommendation.FLAG
            return

        if (self.protein_taxonomy_screen == CommecRecommendation.FLAG or
            self.nucleotide_taxonomy_screen == CommecRecommendation.FLAG):
            self.commec_recommendation = CommecRecommendation.FLAG
            return

        if self.biorisk_screen == CommecRecommendation.WARN:
            self.commec_recommendation = CommecRecommendation.WARN
            return

        if (self.protein_taxonomy_screen == CommecRecommendation.WARN or
            self.nucleotide_taxonomy_screen == CommecRecommendation.WARN):
            self.commec_recommendation = CommecRecommendation.WARN
            return
        # At the moment we only get here when biorisk screen is just run.
        # this will set the global recommend to pass or error, after biorisk is done.
        # We will need earlier checks to maintain the error status of future steps however.
        if self.commec_recommendation == CommecRecommendation.NULL:
            self.commec_recommendation = self.biorisk_screen

@dataclass
class QueryResult:
    '''Container to hold data related to the Query
      used in a commec database screen'''
    query : str = ""
    length : int = 0
    sequence : str = ""
    recommendation : CommecRecommendationContainer = field(default_factory=CommecRecommendationContainer)
    # TODO: Remove summary statistics for first PR?
    summary_info : CommecSummaryStatistics = field(default_factory=CommecSummaryStatistics)
    hits : list[HitDescription] = field(default_factory=list)
    # non_coding_regions : list[list[int]] = field(default_factory=list)

    def get_hit(self, match_name : str) -> HitDescription:
        """
        Searches to see if a match already exists, and returns it, so that it can be modified.
        """
        for hit in self.hits:
            if hit.name == match_name:
                return hit
        return None
    
    def add_new_hit_information(self, new_hit : HitDescription) -> bool:
        """
        Adds a Hit Description to this query, but only adds if the hit is unique, or has a new range.
        Returns True if the hit was not unique, but added to another hit.
        """
        existing_hit = self.get_hit(new_hit.name)
        if not existing_hit:
            self.hits.append(new_hit)
            return False
    
        is_updated : bool = False
        for new_region in new_hit.ranges:
            is_unique_region = True
            for existing_region in existing_hit.ranges:
                if (new_region.query_start == existing_region.query_start and 
                    new_region.query_end == existing_region.query_end):
                    is_unique_region = False
            if is_unique_region:
                is_updated = True
                existing_hit.ranges.append(new_region)

        if is_updated:
            return True

        return False


    def get_flagged_hits(self) -> List[HitDescription]:
        """
        Calculates and returns the list of hits, for all Warnings or Flags.
        Typically used as the regions to check against for benign screens.
        """
        flagged_and_warnings_data = [
        flagged_hit
        for flagged_hit in self.hits if flagged_hit.recommendation.outcome in
        {CommecRecommendation.WARN, CommecRecommendation.FLAG}
        ]
        return flagged_and_warnings_data

    def update(self):
        """
        Call this before exporting to file.
        Sorts the hits based on E-values,
        TODO: Sorts the ranges based on position/E-value (to be confirmed)
        Updates the commec recommendation based on all hits recommendations.
        """
        self.hits.sort(key = lambda hit: hit.get_e_value())
        self.recommendation.update_commec_recommendation(self.hits)

@dataclass
class CommecRunInformation:
    '''Container dataclass to hold general run information for a commec screen '''
    commec_version : str = "0.1.3" # We can programmatically set this eventually.
    json_output_version : str = JSON_COMMEC_FORMAT_VERSION
    biorisk_database_info : SearchToolVersion = field(default_factory=SearchToolVersion)
    protein_database_info : SearchToolVersion = field(default_factory=SearchToolVersion)
    nucleotide_database_info : SearchToolVersion = field(default_factory=SearchToolVersion)
    benign_protein_database_info : SearchToolVersion = field(default_factory=SearchToolVersion)
    benign_rna_database_info : SearchToolVersion = field(default_factory=SearchToolVersion)
    benign_synbio_database_info : SearchToolVersion = field(default_factory=SearchToolVersion)
    time_taken : str = ""
    date_run : str = ""
    # TODO: add other settings / run parameters.

@dataclass
class ScreenResult:
    '''
    Root dataclass to hold all data related to the screening of an individual query by commec.
    '''
    commec_info : CommecRunInformation = field(default_factory = CommecRunInformation)
    queries : list[QueryResult] = field(default_factory=list)

    def format(self):
        ''' Format this ScreenResult as a json string to pass to a standard out if desired.'''
        return str(asdict(self))

    def get_query(self, query_name : str) -> QueryResult:
        """
        Searches for a query, such that it can be updated or read from.
        """
        search_term = query_name
        # Some tools append _X where X is reading frame. Get rid of it for the moment...
        if search_term[-2] == "_":
            search_term = query_name[:-2]

        for data in self.queries:
            if data.query == search_term:
                return data
        return None

    def update(self):
        """
        Propagate update to all children dataclasses.
        """
        for query in self.queries:
            query.update()

    def regions(self) -> Iterator[Tuple[QueryResult, HitDescription, MatchRange]]:
        """
        Helper function, Iterates through all queries, hits, and regions in the ScreenResult object.
        Yields tuples of (query, hit, region).
        """
        for query in self.queries:
            for hit in query.hits:
                for region in hit.ranges:
                    yield query, hit, region

    def hits(self) -> Iterator[Tuple[QueryResult, HitDescription, MatchRange]]:
        """
        helper function, Iterates through all queries, and hits in the ScreenResult object.
        Yields tuples of (query, hit).
        """
        for query in self.queries:
            for hit in query.hits:
                yield query, hit
