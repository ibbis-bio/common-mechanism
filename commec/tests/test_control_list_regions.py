"""
Additional tests for commec.control_list.region.

Complements the existing test_regions.py with tests for load_region_list_data
internals, edge-case inputs, and the REGION_DATA_LUT global.
"""

import json
import pytest

import commec.control_list.region as region_module
from commec.control_list.region import (
    get_regions_set,
    load_region_list_data,
)
from commec.control_list.containers import Region


@pytest.fixture(autouse=True)
def reset_region_lut():
    """Reset the global region lookup between tests."""
    region_module.REGION_DATA_LUT = {}
    yield
    region_module.REGION_DATA_LUT = {}


@pytest.fixture
def region_file(tmp_path):
    """Create a temporary region_definitions.json with an EU group."""
    data = [
        {
            "name": "European Union",
            "acronym": "EU",
            "regions": ["AT", "BE", "DE", "FR"],
        },
    ]
    path = tmp_path / "region_definitions.json"
    path.write_text(json.dumps(data))
    return path


# ---------------------------------------------------------------------------
# load_region_list_data
# ---------------------------------------------------------------------------

def test_load_region_list_data_populates_lut(region_file):
    load_region_list_data(region_file)
    assert "EU" in region_module.REGION_DATA_LUT
    assert region_module.REGION_DATA_LUT["EU"]["name"] == "European Union"
    assert set(region_module.REGION_DATA_LUT["EU"]["regions"]) == {"AT", "BE", "DE", "FR"}


def test_load_region_list_data_missing_file(tmp_path):
    """A non-existent file should leave REGION_DATA_LUT empty without raising."""
    load_region_list_data(tmp_path / "nonexistent.json")
    assert len(region_module.REGION_DATA_LUT) == 0


def test_load_region_list_data_does_not_overwrite_existing(region_file):
    """Loading the same file twice should not overwrite existing entries."""
    load_region_list_data(region_file)
    original_regions = list(region_module.REGION_DATA_LUT["EU"]["regions"])

    # Create a second file that redefines EU with different members.
    alt_data = [{"name": "EU Alt", "acronym": "EU", "regions": ["IT", "ES"]}]
    alt_file = region_file.parent / "alt_regions.json"
    alt_file.write_text(json.dumps(alt_data))
    load_region_list_data(alt_file)

    # Original entry should be preserved.
    assert region_module.REGION_DATA_LUT["EU"]["regions"] == original_regions


def test_load_region_list_data_multiple_groups(tmp_path):
    data = [
        {"name": "Group A", "acronym": "GA", "regions": ["NZ"]},
        {"name": "Group B", "acronym": "GB", "regions": ["AU"]},
    ]
    path = tmp_path / "regions.json"
    path.write_text(json.dumps(data))
    load_region_list_data(path)
    assert "GA" in region_module.REGION_DATA_LUT
    assert "GB" in region_module.REGION_DATA_LUT


# ---------------------------------------------------------------------------
# get_regions_set — edge cases not covered by test_regions.py
# ---------------------------------------------------------------------------

def test_get_regions_set_none_returns_all():
    """None input should map to the special 'all' sentinel."""
    assert get_regions_set(None) == {"all"}


def test_get_regions_set_all_keyword():
    """The literal string 'all' should pass through as-is."""
    assert get_regions_set("all") == {"all"}


def test_get_regions_set_empty_list():
    """An empty list should produce an empty set."""
    assert get_regions_set([]) == set()


def test_get_regions_set_region_object_direct():
    """A Region object whose acronym is a valid alpha-2 code."""
    result = get_regions_set(Region("New Zealand", "NZ"))
    assert "NZ" in result


def test_get_regions_set_list_of_region_objects():
    result = get_regions_set([Region("NZ", "NZ"), Region("AU", "AU")])
    assert result == {"NZ", "AU"}
