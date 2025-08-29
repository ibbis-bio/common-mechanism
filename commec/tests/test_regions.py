import json
import pytest
import logging

from commec.regulation.region import (
    Region,
    get_regions_set,
    load_region_list_data,
)
from commec.utils.logger import setup_console_logging


@pytest.fixture
def region_file(tmp_path):
    """Create a temporary regions.json file with EU and AG groups."""
    data = [
        {
            "name": "European Union",
            "acronym": "EU",
            "regions": [
                "AT", "BE", "BG", "CY", "CZ", "DE", "DK", "EE",
                "ES", "FI", "FR", "GR", "HR", "HU", "IE", "IT",
                "LT", "LU", "LV", "MT", "NL", "PL", "PT", "RO",
                "SE", "SI", "SK"
            ],
        },
        {
            "name": "Australia Group",
            "acronym": "AG",
            "regions": [
                "AR", "AT", "AU", "BE", "BG", "CA", "CH", "CY",
                "CZ", "DE", "DK", "EE", "ES", "FI", "FR", "GB",
                "GR", "HR", "HU", "IE", "IN", "IS", "IT", "JP",
                "KR", "LT", "LU", "LV", "MT", "MX", "NL", "NO",
                "NZ", "PL", "PT", "RO", "SE", "SI", "SK", "TR",
                "UA", "US"
            ],
        },
    ]
    file_path = tmp_path / "regions.json"
    file_path.write_text(json.dumps(data, indent=2))
    return file_path


@pytest.mark.parametrize("acronym,expected", [
    pytest.param(*case) for case in [
        ("AG", {
            "AU", "US", "GB", "DE", "FR", "CA", "JP", "KR", "NZ",
            "NO", "SE", "TR", "UA", "AR", "AT", "BE", "BG", "CH",
            "CY", "CZ", "DK", "EE", "ES", "FI", "GR", "HR", "HU",
            "IE", "IN", "IS", "IT", "LT", "LU", "LV", "MT", "MX",
            "NL", "PL", "PT", "RO", "SI", "SK"
        }), # Custom regions
        (Region("United Kingdom", "GB"), {"GB"}), # From Region
        ([Region("European Union", "EU"), "NZ"], {
            "AT", "BE", "BG", "CY", "CZ", "DE", "DK", "EE", "ES",
            "FI", "FR", "GR", "HR", "HU", "IE", "IT", "LT", "LU",
            "LV", "MT", "NL", "PL", "PT", "RO", "SE", "SI", "SK", "NZ",
        }), # From Region as region custom, with additional country (NZ)
        ("GB", {"GB"}), # Direct
        ("United kingdom", {"GB"}), # Search
        ("United States", {"US"}), # Indirect search
        ("Murica", set()), # Poor search
        (["United States of America", "US"], {"US"}), # Same inputs to set collapse
        (["United States of America", "US", "GB"], {"US", "GB"}), # Different inputs
    ]
])
def test_get_regions_set(region_file, acronym, expected):
    """
    Simple setup, query, compare setup for get_regions_set().
    """
    setup_console_logging(logging.DEBUG)
    load_region_list_data(region_file)
    result = get_regions_set(acronym)
    assert result == expected
