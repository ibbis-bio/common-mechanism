"""Tests for commec.control_list.list_data — module-level state management."""

import pytest
import pandas as pd

import commec.control_list.list_data as ld
from commec.control_list.containers import (
    ControlList,
    ListMode,
    Region,
)
from commec.control_list.initialisation import import_control_lists


@pytest.fixture(autouse=True)
def clean_state():
    """Ensure module state is clean before and after each test."""
    ld.clear()
    yield
    ld.clear()


# ---------------------------------------------------------------------------
# add_control_list
# ---------------------------------------------------------------------------


def test_add_control_list_success():
    cl = ControlList(
        "TestList",
        "Lista Testis",
        "TL",
        "http://example.com",
        Region("NZ", "NZ"),
        ListMode.COMPLIANCE,
        "EXPORT",
    )
    assert ld.add_control_list(cl) is True
    assert "TL" in ld.CONTROL_LISTS
    assert ld.CONTROL_LISTS["TL"] == cl


def test_add_control_list_idempotent():
    """Adding the exact same list twice succeeds (returns True both times)."""
    cl = ControlList(
        "TestList",
        "Lista Testis",
        "TL",
        "http://example.com",
        Region("NZ", "NZ"),
        ListMode.COMPLIANCE,
        "EXPORT",
    )
    assert ld.add_control_list(cl) is True
    assert ld.add_control_list(cl) is True


def test_add_control_list_collision():
    """Different lists sharing an acronym should collide (return False)."""
    cl1 = ControlList(
        "ListA",
        "Lista Prima",
        "TL",
        "http://a.com",
        Region("NZ", "NZ"),
        ListMode.COMPLIANCE,
        "EXPORT",
    )
    cl2 = ControlList(
        "ListB",
        "Lista Secunda",
        "TL",
        "http://b.com",
        Region("AU", "AU"),
        ListMode.COMPLIANCE,
        "EXPORT",
    )
    assert ld.add_control_list(cl1) is True
    assert ld.add_control_list(cl2) is False
    # Original list should be preserved
    assert ld.CONTROL_LISTS["TL"].name == "ListA"


# ---------------------------------------------------------------------------
# Annotation import via CSV files
# ---------------------------------------------------------------------------


def _make_list_dir(
    tmp_path, dirname="testlist", list_info=None, taxids=None, children=None
):
    """Create a minimal valid control-list directory with CSVs."""
    d = tmp_path / dirname
    d.mkdir(parents=True, exist_ok=True)

    info = list_info or [
        {
            "list_name": "TestList",
            "list_acronym": "TL",
            "list_url": "http://example.com",
            "region_name": "New Zealand",
            "region_code": "NZ",
            "use": "EXPORT",
        }
    ]
    pd.DataFrame(info).to_csv(str(d / "list_info.csv"), index=False)

    tx = taxids or [
        {
            "category": "Viruses",
            "display_name": "Test Virus",
            "tax_id": "12345",
            "list_acronym": "TL",
            "list_item": "A.1. Test Virus like organisms",
            "species": "",
            "genus": "",
        }
    ]
    pd.DataFrame(tx).to_csv(str(d / "controlled_taxids.csv"), index=False)

    ch = children or [
        {
            "child_taxid": "99999",
            "controlled_taxid": "12345",
            "child_name": "Son of a Test Virus",
        }
    ]
    pd.DataFrame(ch).to_csv(str(d / "children_of_controlled_taxids.csv"), index=False)

    return d


def test_import_annotations_valid(tmp_path):
    """Annotations from a valid CSV are loaded into CONTROL_LIST_ANNOTATIONS."""
    d = _make_list_dir(tmp_path)
    ld.CONTROL_LIST_ANNOTATIONS.info()
    assert import_control_lists(str(d)) is True
    assert len(ld.CONTROL_LIST_ANNOTATIONS) == 1
    assert ld.CONTROL_LIST_ANNOTATIONS.iloc[0]["tax_id"] == "12345"


def test_import_annotations_minimal_csv_fills_defaults(tmp_path):
    """A CSV with only essential columns should import, filling optional fields."""
    d = _make_list_dir(
        tmp_path,
        taxids=[
            {"tax_id": "100", "list_acronym": "TL", "display_name": "Minimal Org"},
        ],
    )
    assert import_control_lists(str(d)) is True
    row = ld.CONTROL_LIST_ANNOTATIONS.iloc[0]
    assert row["tax_id"] == "100"
    # Columns absent from the CSV should be filled with empty strings
    assert row["species"] == ""
    assert row["genus"] == ""


def test_import_annotations_accumulates_across_lists(tmp_path):
    """Importing from two separate list directories should accumulate annotations."""
    _make_list_dir(
        tmp_path,
        dirname="list_a",
        list_info=[
            {
                "list_name": "L1",
                "list_acronym": "L1",
                "list_url": "u1",
                "region_name": "NZ",
                "region_code": "NZ",
                "use": "EXPORT",
            }
        ],
        taxids=[{"tax_id": "100", "list_acronym": "L1", "name": "Org A"}],
        children=[{"child_taxid": "1", "controlled_taxid": "100", "child_name": "L11"}],
    )
    _make_list_dir(
        tmp_path,
        dirname="list_b",
        list_info=[
            {
                "list_name": "L2",
                "list_acronym": "L2",
                "list_url": "u2",
                "region_name": "AU",
                "region_code": "AU",
                "use": "EXPORT",
            }
        ],
        taxids=[{"tax_id": "200", "list_acronym": "L2", "name": "Org B"}],
        children=[{"child_taxid": "2", "controlled_taxid": "200", "child_name": "L22"}],
    )
    import_control_lists(str(tmp_path / "list_a"))
    import_control_lists(str(tmp_path / "list_b"))
    assert len(ld.CONTROL_LIST_ANNOTATIONS) == 2


# ---------------------------------------------------------------------------
# add_child_lut_data
# ---------------------------------------------------------------------------


def test_add_child_lut_data_missing_columns():
    df = pd.DataFrame([{"child_taxid": "111"}])
    with pytest.raises(ValueError):
        ld.add_child_lut_data(df)


def test_add_child_lut_data_deduplication():
    """Duplicate child-to-parent mappings should be collapsed."""
    df = pd.DataFrame(
        [
            {"child_taxid": "111", "controlled_taxid": "222", "child_name": "Sandy"},
            {"child_taxid": "111", "controlled_taxid": "222", "child_name": "Sandy"},
            {"child_taxid": "333", "controlled_taxid": "444", "child_name": "Dandy"},
            {"child_taxid": "333", "controlled_taxid": "555", "child_name": "Dandy"},
        ]
    )
    ld.add_child_lut_data(df)
    assert len(ld.ACCESSION_MAP) == 3


# ---------------------------------------------------------------------------
# add_ignored_accession_data
# ---------------------------------------------------------------------------


def test_add_ignored_accession_data():
    df = pd.DataFrame(
        [
            {
                "child_taxid": "111",
                "ignored_taxid": "222",
            }
        ]
    )
    ld.add_ignored_accession_data(df)
    assert len(ld.IGNORED_ACCESSION) == 1


def test_add_ignored_accession_data_missing_columns():
    df = pd.DataFrame([{"child_taxid": "111"}])
    with pytest.raises(ValueError):
        ld.add_ignored_accession_data(df)


def test_add_ignored_accession_data_deduplication():
    """Duplicate ignored accessions should be collapsed."""
    df = pd.DataFrame(
        [
            {"child_taxid": "111", "ignored_taxid": "222"},
            {"child_taxid": "111", "ignored_taxid": "222"},
            {"child_taxid": "333", "ignored_taxid": "444"},
        ]
    )
    ld.add_ignored_accession_data(df)
    assert len(ld.IGNORED_ACCESSION) == 2


# ---------------------------------------------------------------------------
# clear
# ---------------------------------------------------------------------------


def test_clear_all():
    """clear() with no argument should reset every piece of module state."""
    ld.add_control_list(
        ControlList(
            "L1",
            "Lista Prima",
            "L1",
            "url",
            Region("NZ", "NZ"),
            ListMode.COMPLIANCE,
            "EXPORT",
        )
    )
    ld.add_control_list_annotations(
        pd.DataFrame([{"list_acronym": "L1", "tax_id": "100"}])
    )
    ld.add_child_lut_data(
        pd.DataFrame(
            [{"child_taxid": "1", "controlled_taxid": "2", "child_name": "Bleargo"}]
        )
    )
    ld.add_ignored_accession_data(
        pd.DataFrame([{"child_taxid": "1", "ignored_taxid": "2"}])
    )

    assert ld.clear() is True
    assert len(ld.CONTROL_LISTS) == 0
    assert len(ld.CONTROL_LIST_ANNOTATIONS) == 0
    assert len(ld.ACCESSION_MAP) == 0
    assert len(ld.IGNORED_ACCESSION) == 0


def test_clear_specific_list():
    """Clearing a specific list by acronym removes only that list and its annotations."""
    ld.add_control_list(
        ControlList(
            "List1",
            "Lista Prima",
            "L1",
            "url1",
            Region("NZ", "NZ"),
            ListMode.COMPLIANCE,
            "EXPORT",
        )
    )
    ld.add_control_list(
        ControlList(
            "List2",
            "Lista Secunda",
            "L2",
            "url2",
            Region("AU", "AU"),
            ListMode.COMPLIANCE,
            "EXPORT",
        )
    )
    ld.add_control_list_annotations(
        pd.DataFrame(
            [
                {"list_acronym": "L1", "tax_id": "100"},
                {"list_acronym": "L2", "tax_id": "200"},
            ]
        )
    )

    assert ld.clear("L1") is True
    assert "L1" not in ld.CONTROL_LISTS
    assert "L2" in ld.CONTROL_LISTS

    remaining_l1 = ld.CONTROL_LIST_ANNOTATIONS[
        ld.CONTROL_LIST_ANNOTATIONS["list_acronym"] == "L1"
    ]
    assert len(remaining_l1) == 0

    remaining_l2 = ld.CONTROL_LIST_ANNOTATIONS[
        ld.CONTROL_LIST_ANNOTATIONS["list_acronym"] == "L2"
    ]
    assert len(remaining_l2) == 1


def test_clear_nonexistent_list():
    """Clearing a list that was never added should return False."""
    assert ld.clear("DOESNOTEXIST") is False
