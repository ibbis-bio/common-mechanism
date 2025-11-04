
"""
The Region module handles the sanitization, and import, of counties
and regions. Supports all current alpha2 country codes through pycountry
module. However also overrides this behaviour with an imported region_definitions.json
data file, which supports general region nesting, such as EU, for the European Union.

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

from .containers import Region

logger = logging.getLogger(__name__)

REGION_DATA_LUT = {}

def load_region_list_data(input_filepath : str | os.PathLike):
    """
    Load the region to country data based on input filepath.
    Input filepath contains a json list containing mapping of 
    two letter codes and names, to a list of affected regions.
    """

    if not os.path.isfile(input_filepath):
        logger.warning("No additional region definitions found at expected location: %s",
        input_filepath)
        return

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
    The input is checked if it is a ...
    * Region : Region, in which case its acronym is returned in the set.
    * Custom Region : str, i.e. EU, returns the set() of all containing alpha-2 codes
    * Arbitrary String : str, uses pycountry to fuzzy search the region, returns the alpha-2 code.

    The use of an arbitrary string is typically unhelpful, unless
    the search term was specific enough to disambigiously identify a single
    region, or small number of regions (<=3 currently). 
    Otherwise, an error is printed, and no country is chosen.
    This is mainly chosen because United States as a search term
    returns 3 [US, UM, VI] (US minor islands, and Virgin Islands)
    """
    global REGION_DATA_LUT

    search_string = region_info

    if region_info == "all" or region_info is None:
        return set(["all"])

    # Handle Region object
    if isinstance(region_info, Region):
        search_string = region_info.acronym

    search_string = str(search_string).strip()

    # Handle Custom region acronyms
    if REGION_DATA_LUT.get(search_string):
        data = REGION_DATA_LUT[search_string]["regions"]
        return set(data)

    # Try 2 letter code retrieval
    if len(search_string) == 2:
        retrieved_country = pc.countries.get(alpha_2 = search_string)
        if retrieved_country:
            return set([search_string])

    # Try 3 letter code retrieval
    if len(search_string) == 3:
        retrieved_country = pc.countries.get(alpha_3 = search_string)
        if retrieved_country:
            return set([retrieved_country.alpha_2])

    # Search with pycountry for alpha_2 codes.
    try:
        countries = pc.countries.search_fuzzy(search_string)
    except LookupError:
        logger.error("Unrecognised region data: \"%s\"", search_string)
        return set()

    # If the search returned with 3 or less outcomes its probably right.
    # This IF statement can be deleted if this behaviour is undesired.
    if len(countries) > 1 and len(countries) <= 3:
        logger.warning("%i countries were returned for input \"%s\".",
                     len(countries), search_string)
        country_list = ",".join([country.alpha_2 for country in countries])
        logger.warning("Countries from search: %s", country_list)
        logger.warning("Using first country: %s", countries[0].name)
        return set([countries[0].alpha_2])

    # If too many values returned, annoy the user.
    if len(countries) > 1:
        logger.error("%i countries were returned for input \"%s\".",
                     len(countries), search_string)
        country_list = ",".join([country.alpha_2 for country in countries])
        logger.error("Countries from search: %s", country_list)
        logger.error("Please refine your input to disambigiously identify a specific country.")
        return set()

    return set([countries[0].alpha_2])



