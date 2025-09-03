
"""
Responsible for the single source of truth - either through ingestion,
or otherwise, where store the relationship between regions that are not 
countries themselves, but groups of countries.

We require this distinction because of the nature of political and policy making
entities not being a true 1:1 with countries. Furthermore, we want the user to input
countries they are interested in from a regulation standpoint, and then get relevant
commec outputs based on this choice.

Someone from Germany would be interested in the policies made that affect
the European Union for example.
"""

import os
import json
import logging
import pycountry as pc

from commec.regulation.containers import Region

logger = logging.getLogger(__name__)

REGION_DATA_LUT = {}

def load_region_list_data(input_filepath : str | os.PathLike):
    """
    Load the region to country data based on input filepath.
    Input filepath contains a json list containing mapping of 
    two letter codes and names, to a list of affected regions.
    """
    region_array = []
    global REGION_DATA_LUT

    with open(input_filepath, encoding = "utf-8") as f:
        region_array = json.load(f)

    for r_data in region_array:
        name = r_data.get("name")
        acronym = r_data.get("acronym")
        region_codes = r_data.get("regions")
        
        if not REGION_DATA_LUT.get(acronym):
            REGION_DATA_LUT[acronym] = {
                "name" : name,
                "regions" : region_codes
            }


def get_regions_set(region_info : str | list[str] | list[Region]) -> set[str]:
    """
    Takes a single, or list of, region information in the form of 
    * alpha-2 codes (i.e. UK, DE, FR)
    * Country names United Kingdom, Germany, France.
    * Region objects Region(name, acronym)
    Validates each entry, and also checks whether the entry is a valid
    commec supported grouping i.e. EU, European Union.
    and returns a set of alpha-2 codes encapulsating all affected regions.
    """
    if not isinstance(region_info, list):
        return set(_return_country_set_from_unknown(region_info))
    
    # Deal with listed inputs.
    flattened = []
    for a in region_info:
        flattened.extend(_return_country_set_from_unknown(a))
    return set(flattened)

def _return_country_set_from_unknown(region_info : str | Region = "") -> set[str]:
    """
    Given a String or Region object, return a set of alpha-2 country codes.
    The input is checked if it is a 
    * Region, in which case its acronym is returned in the set.
    * Custom Region i.e. EU, returns the set() of all containing alpha-2 codes
    * Arbitrary String, uses pycountry to fuzzy search the region, and returns the alpha-2 code.
    """
    global REGION_DATA_LUT

    search_string = region_info

    # Handle Region object
    if isinstance(region_info, Region):
        search_string = region_info.acronym

    # Handle Custom region acronyms
    if REGION_DATA_LUT.get(search_string):
        data = REGION_DATA_LUT[search_string]["regions"]
        return set(data)

    # Try 2 letter code retrieval
    if len(search_string) == 2:
        retrieved_country = pc.countries.get(alpha_2 = search_string)
        if retrieved_country:
            return set([search_string])

    # Search with pycountry for alpha_2 codes.
    try:
        countries = pc.countries.search_fuzzy(search_string)
    except(LookupError):
        logger.error("Unrecognised region data: \"%s\"", search_string)
        return set()

    if len(countries) > 1:
        logger.debug("%i countries were returned for input \"%s\".",
                     len(countries), search_string)

    return [countries[0].alpha_2]



