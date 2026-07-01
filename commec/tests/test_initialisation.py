"""Tests for commec.control_list.initialisation — data import pipeline."""

import os

import pytest
import pandas as pd

import commec.control_list.list_data as ld
from commec.control_list.containers import (
    Accession,
    CategoryType,
    ControlList,
    ListMode,
    Region,
)
from commec.control_list.initialisation import (
    import_control_lists,
    tidy_control_list_data,
    update_regional_context,
)


# ---------------------------------------------------------------------------
# Fixtures & Helpers
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def clean_state():
    ld.clear()
    yield
    ld.clear()


def _write_csv(directory, filename, rows):
    """Write a list of dicts to a CSV inside *directory*."""
    pd.DataFrame(rows).to_csv(os.path.join(str(directory), filename), index=False)


def _make_valid_list_dir(tmp_path, dirname="testlist", extra_lists=None,
                         extra_taxids=None, extra_children=None,
                         ignored=None, taxids=None, list_info=None):
    """
    Create a minimal valid control-list directory containing:
      list_info.csv, controlled_taxids.csv, children_of_controlled_taxids.csv
    and optionally ignored_taxids.csv.
    Returns the directory Path.
    """
    d = tmp_path / dirname
    d.mkdir(parents=True, exist_ok=True)

    info = list_info or [{
        "list_name": "TestList", "list_acronym": "TL",
        "list_url": "http://example.com",
        "region_name": "New Zealand", "region_code": "NZ", "use": "EXPORT",
    }]
    if extra_lists:
        info.extend(extra_lists)
    _write_csv(d, "list_info.csv", info)

    tx = taxids or [{
        "category": "Viruses", "display_name": "Test Virus",
        "tax_id": "11320", "list_acronym": "TL",
        "species": "TestSpecies",
        "genus": "TestGenus",
    }]
    if extra_taxids:
        tx.extend(extra_taxids)
    _write_csv(d, "controlled_taxids.csv", tx)

    children = extra_children or [
        {"child_taxid": "99999", "controlled_taxid": "11320", "child_name" : "Test Viruses first born son"},
    ]
    _write_csv(d, "children_of_controlled_taxids.csv", children)

    if ignored is not None:
        _write_csv(d, "ignored_taxids.csv", ignored)

    return d


# ---------------------------------------------------------------------------
# import_control_lists — required / optional files
# ---------------------------------------------------------------------------

def test_import_control_lists_valid_directory(tmp_path):
    d = _make_valid_list_dir(tmp_path)
    assert import_control_lists(str(d)) is True
    assert "TL" in ld.CONTROL_LISTS
    assert len(ld.CONTROL_LIST_ANNOTATIONS) >= 1
    assert len(ld.ACCESSION_MAP) >= 1


def test_import_control_lists_missing_controlled_taxids(tmp_path):
    d = tmp_path / "incomplete"
    d.mkdir()
    _write_csv(d, "list_info.csv", [{
        "list_name": "TL", "list_acronym": "TL", "list_url": "url",
        "region_name": "NZ", "region_code": "NZ", "use": "EXPORT",
    }])
    _write_csv(d, "children_of_controlled_taxids.csv", [
        {"child_taxid": "1", "controlled_taxid": "2"},
    ])
    # controlled_taxids.csv is missing
    assert import_control_lists(str(d)) is False


def test_import_control_lists_missing_children_csv(tmp_path):
    d = tmp_path / "nochildren"
    d.mkdir()
    _write_csv(d, "list_info.csv", [{
        "list_name": "TL", "list_acronym": "TL", "list_url": "url",
        "region_name": "NZ", "region_code": "NZ", "use": "EXPORT",
    }])
    _write_csv(d, "controlled_taxids.csv", [
        {"tax_id": "100", "list_acronym": "TL", "display_name": "Org"},
    ])
    # children_of_controlled_taxids.csv is missing
    assert import_control_lists(str(d)) is False


def test_import_control_lists_nonexistent_directory():
    assert import_control_lists("/this/path/does/not/exist") is False


def test_import_control_lists_optional_ignored_accessions(tmp_path):
    """Import should succeed when ignored_taxids.csv is absent."""
    d = _make_valid_list_dir(tmp_path)
    assert not os.path.isfile(os.path.join(str(d), "ignored_taxids.csv"))
    assert import_control_lists(str(d)) is True


def test_import_control_lists_with_ignored_accessions(tmp_path):
    d = _make_valid_list_dir(tmp_path,
                             ignored=[{"child_taxid": "50", "ignored_taxid": "100"}])
    assert import_control_lists(str(d)) is True
    assert len(ld.IGNORED_ACCESSION) == 1


# ---------------------------------------------------------------------------
# Comma-separated values in CSVs
# ---------------------------------------------------------------------------

def test_import_annotations_comma_separated_taxids(tmp_path):
    """Comma-separated tax_id values should be exploded into separate rows."""
    d = _make_valid_list_dir(tmp_path, taxids=[
        {"tax_id": "100,200,300", "list_acronym": "TL", "display_name": "Multi"},
    ])
    assert import_control_lists(str(d)) is True
    assert len(ld.CONTROL_LIST_ANNOTATIONS) == 3


def test_import_annotations_comma_separated_list_acronyms(tmp_path):
    """Comma-separated list_acronym values should be exploded into separate rows."""
    d = _make_valid_list_dir(
        tmp_path,
        list_info=[
            {"list_name": "L1", "list_acronym": "L1", "list_url": "url1",
             "region_name": "NZ", "region_code": "NZ", "use": "EXPORT"},
            {"list_name": "L2", "list_acronym": "L2", "list_url": "url2",
             "region_name": "AU", "region_code": "AU", "use": "EXPORT"},
        ],
        taxids=[
            {"tax_id": "100", "list_acronym": "L1, L2", "display_name": "Shared Organism"},
        ],
    )
    assert import_control_lists(str(d)) is True
    assert len(ld.CONTROL_LIST_ANNOTATIONS) == 2


def test_import_annotations_invalid_acronym_filtered(tmp_path):
    """Annotations referencing unknown list acronyms should be silently dropped."""
    d = _make_valid_list_dir(tmp_path, taxids=[
        {"tax_id": "100", "list_acronym": "TL", "display_name": "Valid"},
        {"tax_id": "200", "list_acronym": "NOPE", "display_name": "Invalid"},
    ])
    assert import_control_lists(str(d)) is True
    # Only the row for "TL" should survive.
    assert len(ld.CONTROL_LIST_ANNOTATIONS) == 1


# ---------------------------------------------------------------------------
# update_regional_context
# ---------------------------------------------------------------------------

def test_update_regional_context_matching_region():
    ld.add_control_list(ControlList(
        "L1", "Lista Prima", "L1", "url", Region("New Zealand", "NZ"),
        ListMode.CONDITIONAL_NUM, "EXPORT"))
    update_regional_context({"NZ"})
    assert ld.CONTROL_LISTS["L1"].status == ListMode.COMPLIANCE


def test_update_regional_context_non_matching():
    ld.add_control_list(ControlList(
        "L1", "Lista Prima", "L1", "url", Region("New Zealand", "NZ"),
        ListMode.COMPLIANCE, "EXPORT"))
    update_regional_context({"AU"})
    assert ld.CONTROL_LISTS["L1"].status == ListMode.CONDITIONAL_NUM


def test_update_regional_context_all_keyword():
    """The 'all' sentinel should force every list to COMPLIANCE."""
    ld.add_control_list(ControlList(
        "L1", "Lista Prima", "L1", "url1", Region("New Zealand", "NZ"),
        ListMode.CONDITIONAL_NUM, "EXPORT"))
    ld.add_control_list(ControlList(
        "L2", "Lista Secunda", "L2", "url2", Region("Australia", "AU"),
        ListMode.CONDITIONAL_NUM, "EXPORT"))
    update_regional_context({"all"})
    assert ld.CONTROL_LISTS["L1"].status == ListMode.COMPLIANCE
    assert ld.CONTROL_LISTS["L2"].status == ListMode.COMPLIANCE


def test_update_regional_context_none_defaults_to_all():
    """Passing None should behave like 'all' (no restriction)."""
    ld.add_control_list(ControlList(
        "L1", "Lista Prima", "L1", "url", Region("New Zealand", "NZ"),
        ListMode.CONDITIONAL_NUM, "EXPORT"))
    update_regional_context(None)
    assert ld.CONTROL_LISTS["L1"].status == ListMode.COMPLIANCE


def test_update_regional_context_custom_alternative_mode():
    ld.add_control_list(ControlList(
        "L1", "Lista Prima", "L1", "url", Region("New Zealand", "NZ"),
        ListMode.COMPLIANCE, "EXPORT"))
    update_regional_context({"AU"}, alternative_mode=ListMode.IGNORE)
    assert ld.CONTROL_LISTS["L1"].status == ListMode.IGNORE


# ---------------------------------------------------------------------------
# tidy_control_list_data
# ---------------------------------------------------------------------------

def test_tidy_creates_accession_index():
    ld.add_control_list(ControlList(
        "L1", "Lista Prima", "L1", "url", Region("New Zealand", "NZ"), ListMode.COMPLIANCE, "EXPORT"))
    ld.add_control_list_annotations(pd.DataFrame([{
        "list_acronym": "L1", "tax_id": "12345",
        "display_name": "Virus", "category": "Viruses", "lineage": "",
    }]))
    tidy_control_list_data()
    assert isinstance(ld.CONTROL_LIST_ANNOTATIONS.index[0], Accession)


def test_tidy_removes_invalid_accessions():
    """Entries with empty tax_id should be dropped during tidying."""
    ld.add_control_list(ControlList(
        "L1", "Lista Prima", "L1", "url", Region("New Zealand", "NZ"), ListMode.COMPLIANCE, "EXPORT"))
    ld.add_control_list_annotations(pd.DataFrame([
        {"list_acronym": "L1", "tax_id": "12345",
         "display_name": "Valid", "category": "Viruses"},
        {"list_acronym": "L1", "tax_id": "",
         "display_name": "BadEntry", "category": "Viruses"},
    ]))
    print("Current Control List Data before tidying: %s", ld.CONTROL_LIST_ANNOTATIONS)
    tidy_control_list_data()
    print("Control List Data after tidying: %s", ld.CONTROL_LIST_ANNOTATIONS)
    assert len(ld.CONTROL_LIST_ANNOTATIONS) == 1


def test_tidy_deduplicates_by_list_and_accession():
    """Exact duplicates on (list_acronym, accession) should be collapsed."""
    ld.add_control_list(ControlList(
        "L1", "Lista Prima", "L1", "url", Region("New Zealand", "NZ"), ListMode.COMPLIANCE, "EXPORT"))
    ld.add_control_list_annotations(pd.DataFrame([
        {"list_acronym": "L1", "tax_id": "100",
         "display_name": "Same", "category": "Viruses"},
        {"list_acronym": "L1", "tax_id": "100",
         "display_name": "Same", "category": "Viruses"},
    ]))
    tidy_control_list_data()
    assert len(ld.CONTROL_LIST_ANNOTATIONS) == 1


def test_tidy_invalid_category_defaults_to_none():
    ld.add_control_list(ControlList(
        "L1", "Lista Prima", "L1", "url", Region("New Zealand", "NZ"), ListMode.COMPLIANCE, "EXPORT"))
    ld.add_control_list_annotations(pd.DataFrame([{
        "list_acronym": "L1", "tax_id": "100",
        "display_name": "Test", "category": "InvalidCategory",
    }]))
    tidy_control_list_data()
    assert ld.CONTROL_LIST_ANNOTATIONS.iloc[0]["category"] == CategoryType.NONE


def test_tidy_valid_category_preserved():
    ld.add_control_list(ControlList(
        "L1", "Lista Prima", "L1", "url", Region("New Zealand", "NZ"), ListMode.COMPLIANCE, "EXPORT"))
    ld.add_control_list_annotations(pd.DataFrame([{
        "list_acronym": "L1", "tax_id": "100",
        "display_name": "Test", "category": "Bacteria",
    }]))
    tidy_control_list_data()
    assert ld.CONTROL_LIST_ANNOTATIONS.iloc[0]["category"] == CategoryType.BACTERIA


# ---------------------------------------------------------------------------
# _import_control_list_info — invalid 'use' column
# ---------------------------------------------------------------------------

def test_import_list_info_invalid_use_column(tmp_path):
    """An invalid 'use' value should prevent the list from being imported."""
    d = _make_valid_list_dir(
        tmp_path, dirname="bad_use",
        list_info=[{
            "list_name": "Bad", "list_acronym": "BAD",
            "list_url": "url", "region_name": "NZ",
            "region_code": "NZ", "use": "INVALID_USE",
        }],
    )
    import_control_lists(str(d))
    assert "BAD" not in ld.CONTROL_LISTS
