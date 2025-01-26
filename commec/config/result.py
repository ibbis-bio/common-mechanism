#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
'''
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
        [QueryResults]
            recommendation
            [HitResults]:
                recommendation
'''

from dataclasses import dataclass, asdict, field
from typing import List, Iterator, Tuple
from enum import StrEnum
from importlib.metadata import version, PackageNotFoundError
from commec.tools.search_handler import SearchToolVersion
from commec.config.query import Query

try:
    COMMEC_VERSION = version("commec")
except PackageNotFoundError:
    COMMEC_VERSION = "error"
# Seperate versioning for the output JSON.
JSON_COMMEC_FORMAT_VERSION = "0.1"

class Recommendation(StrEnum):
    """
    All possible outputs from commec for a query screen, or individual hit.
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
            Recommendation.NULL: 0,
            Recommendation.SKIP: 1,
            Recommendation.PASS: 2,
            Recommendation.CLEARED_WARN: 3,
            Recommendation.CLEARED_FLAG: 4,
            Recommendation.WARN: 5,
            Recommendation.FLAG: 6,
            Recommendation.ERROR: 7,
        }
        return order[self]

    def clear(self):
        """ Convert a WARN or FLAG into its cleared counterpart, and return that. """
        if self == Recommendation.WARN:
            return Recommendation.CLEARED_WARN
        if self == Recommendation.FLAG:
            return Recommendation.CLEARED_FLAG
        return self

def compare(a : Recommendation, b : Recommendation):
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
    BIORISK = 'Biorisk Screen'
    TAXONOMY_NT = "Nucleotide Taxonomy Screen"
    TAXONOMY_AA = "Protein Taxonomy Screen"
    BENIGN_PROTEIN = "Benign Protein Screen"
    BENIGN_RNA = "Benign RNA Screen"
    BENIGN_SYNBIO = "Benign SynBio Screen"

@dataclass
class HitRecommendationContainer:
    """ 
    Summarises the recommendation and step for a single Hit.
    Pairs a recommendation with a Screening step.
    """
    outcome : Recommendation = Recommendation.NULL
    from_step : ScreenStep = field(default_factory=ScreenStep)

@dataclass
class MatchRange:
    """
    Container for coordinate information of where hits match to a query.
    """
    e_value : float = 0.0
    # percent identity?
    match_start : int = 0
    match_end : int = 0
    query_start : int = 0
    query_end : int = 0

    # TODO: Add frame, as QueryStart and QueryEnd should be in Frame0 NT coords.

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
class HitResult:
    """ 
    Container for all information regarding a single hit with a range(s), to a single query.
    A Hit is any outcome from any step during Commec Screening, which can be mapped onto the Query.
    """
    recommendation : HitRecommendationContainer = field(default_factory=HitRecommendationContainer)
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
class QueryRecommendationContainer:
    """
    Summarises the recommendations across all hits for a single query.
    """
    commec_recommendation : Recommendation = Recommendation.NULL
    biorisk_screen : Recommendation = Recommendation.NULL
    protein_taxonomy_screen : Recommendation = Recommendation.NULL
    nucleotide_taxonomy_screen : Recommendation = Recommendation.NULL
    benign_screen : Recommendation = Recommendation.NULL

    def update_commec_recommendation(self):
        """
        Parses the current state of the json,
        and updates the global commec recommendation.

        TODO: Consider the value of iterating over all HitResults (provided as input)
        and updating things with stricter logic.
        """
        if self.benign_screen in {Recommendation.CLEARED_FLAG,
                                  Recommendation.CLEARED_WARN,
                                  Recommendation.PASS}:
            self.commec_recommendation = self.benign_screen
            return

        if self.biorisk_screen == Recommendation.FLAG:
            self.commec_recommendation = Recommendation.FLAG
            return

        if (self.protein_taxonomy_screen == Recommendation.FLAG or
            self.nucleotide_taxonomy_screen == Recommendation.FLAG):
            self.commec_recommendation = Recommendation.FLAG
            return

        if self.biorisk_screen == Recommendation.WARN:
            self.commec_recommendation = Recommendation.WARN
            return

        if (self.protein_taxonomy_screen == Recommendation.WARN or
            self.nucleotide_taxonomy_screen == Recommendation.WARN):
            self.commec_recommendation = Recommendation.WARN
            return
        # At the moment we only get here when biorisk screen is just run.
        # this will set the global recommend to pass or error, after biorisk is done.
        # We will need earlier checks to maintain the error status of future steps however.
        if self.commec_recommendation == Recommendation.NULL:
            self.commec_recommendation = self.biorisk_screen


@dataclass
class QueryResult:
    '''
    Container to hold screening result data pertinant to a single Query
    '''
    query : str = ""
    length : int = 0
    sequence : str = ""
    recommendation : QueryRecommendationContainer = field(default_factory=QueryRecommendationContainer)
    hits : dict[str, HitResult] = field(default_factory=dict)

    def get_hit(self, match_name : str) -> HitResult:
        """
        Searches to see if a match already exists, and returns it, so that it can be modified.
        """
        return self.hits.get(match_name)
        for hit in self.hits:
            if hit.name == match_name:
                return hit
        return None
    
    def add_new_hit_information(self, new_hit : HitResult) -> bool:
        """
        Adds a Hit Description to this query, but only adds if the hit is unique, or has a new range.
        Returns True if the hit was not unique, but added unique info to the hit.
        """
        existing_hit = self.hits.get(new_hit.name)
        if not existing_hit:
            self.hits[new_hit.name] = new_hit
            return False
    
        hits_is_updated : bool = False
        for new_region in new_hit.ranges:
            is_unique_region = True
            for existing_region in existing_hit.ranges:
                if (new_region.query_start == existing_region.query_start and 
                    new_region.query_end == existing_region.query_end):
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
        for flagged_hit in self.hits.values() if flagged_hit.recommendation.outcome in
        {Recommendation.WARN, Recommendation.FLAG}
        ]
        return flagged_and_warnings_data

    def update(self):
        """
        Call this before exporting to file.
        Sorts the hits based on E-values,
        TODO: Sorts the ranges based on position/E-value (to be confirmed)
        Updates the commec recommendation based on all hits recommendations.
        """
        # A rare instance where we want our dictionary to be sorted - 
        # So we convert it to a list, and recreate a dictionary - which maintains the order.
        # Thankfully we only do this once per Query at the end.
        #sorted_items_desc = sorted(self.hits.items(), key=lambda item: item.get_e_value(), reverse=True)
        #self.hits = dict(sorted_items_desc)
        self.recommendation.update_commec_recommendation()

@dataclass
class ScreenRunInfo:
    '''Container dataclass to hold general run information for a commec screen '''
    commec_version : str = str(COMMEC_VERSION)
    json_output_version : str = JSON_COMMEC_FORMAT_VERSION
    biorisk_database_info : SearchToolVersion = field(default_factory=SearchToolVersion)
    protein_database_info : SearchToolVersion = field(default_factory=SearchToolVersion)
    nucleotide_database_info : SearchToolVersion = field(default_factory=SearchToolVersion)
    benign_protein_database_info : SearchToolVersion = field(default_factory=SearchToolVersion)
    benign_rna_database_info : SearchToolVersion = field(default_factory=SearchToolVersion)
    benign_synbio_database_info : SearchToolVersion = field(default_factory=SearchToolVersion)
    time_taken : str = ""
    date_run : str = ""

@dataclass
class ScreenResult:
    '''
    Root dataclass to hold all data related to the screening of an individual query by commec.
    '''
    commec_info : ScreenRunInfo = field(default_factory = ScreenRunInfo)
    queries : dict[str, QueryResult] = field(default_factory=dict)

    def format(self):
        ''' Format this ScreenResult as a json string to pass to a standard out if desired.'''
        return str(asdict(self))

    def get_query(self, query_name : str) -> QueryResult:
        """
        Searches for a query, such that it can be updated or read from.
        """
        search_term = Query.create_id(query_name)
        return self.queries.get(search_term)
    
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
