#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
'''
   Set of tools for retrieving and storing information important to screen
    outputs. Information is stored as a structure of dataclasses, and are 
    converted between the dataclass / dict / json_file as required. The 
    conversions are done dynamically, and it is recommended to only use and 
    interact with the dataclasses only, to maintain version format, and not 
    create erroneous outputs to the JSON which wont be read back in. This
    ensures an expected i/o behaviour.

   In this way, the JSON object serves as a common state, that can be updated
    whilst not being temporally appended like a log file i.e. .screen file.

    Example JSON structure:
    {
        commec_run_info : {
            json / commec / database versions etc.
            time run, date run etc.
        }
    }
'''

# Consider whether this can get away with being part of config. rename to IO config?

import json
import string
import os
from dataclasses import dataclass, asdict, fields, field, is_dataclass
from typing import Dict, Type, get_origin, Any, get_args
from enum import StrEnum

# Seperate versioning for the output JSON.
JSON_COMMEC_FORMAT_VERSION = "0.1"

class CommecRecomendation(StrEnum):
    """
    All possible recommendation outputs from commec for a query.
    """
    NULL = '-' # This was not set.
    SKIP = 'Skip' # Intentionally skipped this step.
    PASS = 'Pass' # Commec has approved this query at this step.
    WARN = 'Warn' # This query may be suspicious...
    FLAG = 'Flag' # Commec has flagged this query.
    ERROR = 'Error' # An error occured, such that this step failed to run.

    # I don't like the below, should we just use FLAG and PASS perhaps.
    TRUE = 'true' # Useful for organism identification outputs.
    FALSE = 'false' # Useful for organism identification outputs.

@dataclass
class CommecRecomendationContainer():
    """
    Summarises the recommendations at each step of the commec screen process.
    """
    commec_recommendation : CommecRecomendation = CommecRecomendation.NULL
    biorisk_screen : CommecRecomendation = CommecRecomendation.NULL
    protein_taxonomy_screen : CommecRecomendation = CommecRecomendation.NULL
    nucleotide_taxonomy_screen : CommecRecomendation = CommecRecomendation.NULL
    virus_flag : CommecRecomendation = CommecRecomendation.NULL
    bacteria_flag : CommecRecomendation = CommecRecomendation.NULL
    eukaryote_flag : CommecRecomendation = CommecRecomendation.NULL
    benign_screen : CommecRecomendation = CommecRecomendation.NULL

    def update_commec_recommendation(self):
        """ 
        Parses the current state of the json, 
        and updates the global commec recommendation.
        """

        if self.biorisk_screen == CommecRecomendation.FLAG:
            self.commec_recommendation = CommecRecomendation.FLAG
            return

        if (self.protein_taxonomy_screen == CommecRecomendation.FLAG
            and self.benign_screen == CommecRecomendation.FLAG):
            self.commec_recommendation = CommecRecomendation.FLAG
            return
        
        if (self.nucleotide_taxonomy_screen == CommecRecomendation.FLAG
            and self.benign_screen == CommecRecomendation.FLAG):
            self.commec_recommendation = CommecRecomendation.FLAG
            return
        
        if self.biorisk_screen == CommecRecomendation.WARN:
            self.commec_recommendation = CommecRecomendation.WARN
            return

        if self.commec_recommendation == CommecRecomendation.NULL:
            self.commec_recommendation = self.biorisk_screen


@dataclass
class CommecSummaryStatistics:
    """
    Summary of some useful statistics for rapid intepretation of complicated screen outputs.
    """
    biorisk_toxin_percentage : int = 0
    percentage_of_regulated_for_taxonomy_best_matches : int = 0
    human_pathogen_flag_percentage : int = 0

@dataclass
class MatchRange:
    '''Container for information on how a match maps to a query, match direction
      is +1, or positive, if the query and match are in the same direction, and 
      -1 if they are inverse, and 0 if not set.'''
    match_start : int = 0
    match_end : int = 0
    query_start : int = 0
    query_end : int = 0
    query_direction : int = 0

@dataclass
class BenignData:
    '''Container to hold data related to the benign status of a taxonomy search'''
    benign : bool = False
    benign_match_type : str = "" # CMscan, Blast, HMM
    percent_identity_e_value : float = 0.0

@dataclass
class MatchFields:
    '''Container to hold information for a match to the query identified in a commec screen. 
    molecule_alphabet is expected to contain either 'aa' for aminoacids, or 'nt' for nucelotides.
    range contains information pertaining to over what range the match occured.'''
    description : str = ""
    molecule_alphabet : str = ""
    id : int = 0
    taxon : str = ""
    kingdom : str = ""
    ranges : list[MatchRange] = field(default_factory = list)
    is_benign : BenignData = field(default_factory = BenignData)
    # etc etc

@dataclass
class BioRisk:
    '''Container to hold information for a match to the query identified as a potential biorisk'''
    target_name : str = ""
    description : str = ""
    regulated : bool = True # This may already be encoded in the list that the instance exists in.
    regulated_info : str = ""
    range : list[MatchRange] = field(default_factory = list)
    # etc etc

@dataclass
class BioRiskData:
    '''
    Container dataclass for a list of matches to biorisks 
    identified from a commec database screen.
    '''
    regulated_genes : list[BioRisk] = field(default_factory = list)
    virulance_factors : list[BioRisk] = field(default_factory = list)

    def get_existing(self, match_description : str):
        """
        Searches to see if a match already exists, and returns it, so that it can be modified.
        """
        for regulated_match in self.regulated_genes:
            if regulated_match.description == match_description:
                return regulated_match
        for virulance_factor in self.virulance_factors:
            if virulance_factor.description == match_description:
                return virulance_factor
        return None

    def get_existing_regulated_gene(self, match_description : str):
        """
        Searches to see if a match already exists, and returns it, so that it can be modified.
        """
        for regulated_match in self.regulated_genes:
            if regulated_match.description == match_description:
                return regulated_match
        return None

    def get_existing_virulance_factor(self, match_description : str):
        """
        Searches to see if a match already exists, and returns it, so that it can be modified.
        """
        for virulance_factor in self.virulance_factors:
            if virulance_factor.description == match_description:
                return virulance_factor
        return None

@dataclass
class TaxonomyData:
    '''
    Container dataclass for a list of matches of taxonomy hits, 
    identified from a commec database screen.
    '''
    is_regulated : bool
    regulation_agency : str = ""
    matches : list[MatchFields] = field(default_factory = list)

    def get_existing_match(self, match_description : str):
        """
        Searches to see if a match already exists, and returns it, so that it can be modified.
        """
        for match in self.matches:
            if match.description == match_description:
                return match
        return None

@dataclass
class QueryData:
    '''Container to hold data related to the Query
      used in a commec database screen'''
    name : str = ""
    length : int = 0
    sequence : str = ""
    recommendation : CommecRecomendationContainer = field(default_factory=CommecRecomendationContainer)
    summary_info : CommecSummaryStatistics = field(default_factory=CommecSummaryStatistics)
    biorisks : BioRiskData = field(default_factory = BioRiskData)
    taxonomies : list[TaxonomyData] = field(default_factory = list)

@dataclass
class CommecRunInformation:
    '''Container dataclass to hold general run information for a commec screen '''
    commec_version : str = "0.1.3" # We can programmatically set this eventually.
    json_output_version : str = JSON_COMMEC_FORMAT_VERSION
    biorisk_database_info : str = ""
    protein_database_info : str = ""
    nucleotide_database_info : str = ""
    benign_database_info : str = ""
    time_taken : str = ""
    date_run : str = ""
    # add other settings / run parameters.

@dataclass
class ScreenData:
    ''' Root dataclass to hold all data related to the screening of an individual query by commec.'''
    commec_info : CommecRunInformation = field(default_factory = CommecRunInformation)
    queries : list[QueryData] = field(default_factory = list)

    def format(self):
        ''' Format this ScreenData as a json string to pass to a standard out if desired.'''
        return str(asdict(self))
    
    def get_query(self, query_name : str) -> QueryData:
        """
        Searches for a query, such that it can be updated or read from.
        """
        search_term = query_name
        # Some tools append _X where X is reading frame. Get rid of it for the moment...
        if search_term[-2] == "_":
            search_term = query_name[:-2]

        for data in self.queries:
            if data.name == search_term:
                return data
        return None


# The above could be moved to a custom .py script for variable importing under version control.

def encode_screen_data_to_json(input_screendata: ScreenData, output_json_filepath: string = "output.json") -> None:
    ''' Converts a ScreenData class object into a JSON file at the given filepath.'''
    with open(output_json_filepath, "w", encoding="utf-8") as json_file:
        json.dump(asdict(input_screendata), json_file, indent=4)

def encode_dict_to_screen_data(input_dict : dict) -> ScreenData:
    ''' Converts a dictionary into a ScreenData object, 
    any keys within the dictionary not part of the ScreenData format are lost.
    any missing information will be simple set as defaults.'''
    return dict_to_dataclass(ScreenData, input_dict)



# Convert the dictionary back to the dataclass or list of dataclass
def dict_to_dataclass(cls: Type, data: Dict[str, Any]) -> Any:
    ''' 
    Convert a dict, into appropriate dataclass, or list of dataclass, 
    invalid keys to the dataclass structure are ignored.
    '''
    # Prepare a dictionary for filtered data
    filtered_data = {}

    for f in fields(cls):
        field_name = f.name
        field_type = f.type
        if field_name in data:
            field_value = data[field_name]

            # Check if the field is a dataclass
            if is_dataclass(field_type):
                filtered_data[field_name] = dict_to_dataclass(field_type, field_value)
                continue

            # Check if the field is a list of dataclasses
            if get_origin(field_type) is list:
                item_type = get_args(field_type)[0]
                if is_dataclass(item_type) and isinstance(field_value, list):
                    filtered_data[field_name] = [
                        dict_to_dataclass(item_type, item) for item in field_value
                        if isinstance(item, dict) 
                            and any(key in {f.name for f in fields(item_type)} for key in item.keys())
                            or isinstance(item, item_type)]
                    continue
                filtered_data[field_name] = field_value
                continue

            # Handle other field types
            filtered_data[field_name] = field_value

    # Create an instance of the dataclass with the filtered data
    return cls(**filtered_data)

def get_screen_data_from_json(input_json_filepath: string) -> ScreenData:
    ''' Loads a JSON file from given filepath and returns 
    a populated ScreenData object from its contents. If the file does not
    exist, then returns a new screen data object.'''
    if not os.path.exists(input_json_filepath):
        return ScreenData()

    json_string : str
    with open(input_json_filepath, "r", encoding="utf-8") as json_file:
        # Read the file contents as a string
        json_string = json_file.read()
    my_data : dict = json.loads(json_string)

    # Check version of imported json.
    input_version = my_data["commec_info"]["json_output_version"]
    if not input_version == JSON_COMMEC_FORMAT_VERSION:
        raise RuntimeError("""Version difference between input (v.{input_version}) and
                            expected (v.{JSON_COMMEC_FORMAT_VERSION}) JSON state file: {input_json_filepath}""")
    return encode_dict_to_screen_data(my_data)