
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
import pycountry as pc

_region_to_country_data : dict[str,list[str]] = {
    "EU": ["DE", "FR"],
}

def load_region_list_data(input_filepath : str | os.PathLike):
    """
    Load the region to country data from a single source of truth...
    To be decided, source is hardcoded for now.
    """
    ...