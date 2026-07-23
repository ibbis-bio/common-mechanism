"""Tests for commec.control_list.control_list — main public API functions."""

import json
import os

import pytest
import pandas as pd

import commec.control_list.list_data as ld
import commec.control_list.region as region_module
from commec.control_list.containers import (
    ControlList,
    ListMode,
    Region,
)
from commec.control_list.control_list import (
    get_control_lists,
    get_regulation,
    import_data,
    is_regulated,
)
from commec.control_list.initialisation import tidy_control_list_data


# ---------------------------------------------------------------------------
# Fixtures & Helpers
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def clean_state():
    ld.clear()
    region_module.REGION_DATA_LUT = {}
    yield
    ld.clear()
    region_module.REGION_DATA_LUT = {}


def _build_db(base_path, lists):
    """
    Build a control-list database directory for ``import_data``.

    *base_path* : root directory (will contain ``region_definitions.json``)
    *lists*     : list of dicts, each with keys:
        folder_name, list_info (list[dict]), taxids, children,
        and optionally ignored.
    """
    os.makedirs(str(base_path), exist_ok=True)
    region_defs = os.path.join(str(base_path), "region_definitions.json")
    if not os.path.isfile(region_defs):
        with open(region_defs, "w") as f:
            json.dump([], f)

    for entry in lists:
        folder = os.path.join(str(base_path), entry["folder_name"])
        os.makedirs(folder, exist_ok=True)
        pd.DataFrame(entry["list_info"]).to_csv(
            os.path.join(folder, "list_info.csv"), index=False
        )
        pd.DataFrame(entry.get("taxids", [])).to_csv(
            os.path.join(folder, "controlled_taxids.csv"), index=False
        )
        pd.DataFrame(entry.get("children", [])).to_csv(
            os.path.join(folder, "children_of_controlled_taxids.csv"), index=False
        )
        if "ignored" in entry:
            pd.DataFrame(entry["ignored"]).to_csv(
                os.path.join(folder, "ignored_taxids.csv"), index=False
            )


_DEFAULT_LIST = {
    "folder_name": "mylist",
    "list_info": [
        {
            "list_name": "TestList",
            "list_acronym": "TL",
            "list_url": "http://example.com",
            "region_name": "New Zealand",
            "region_code": "NZ",
            "use": "EXPORT",
        }
    ],
    "taxids": [
        {
            "tax_id": "100",
            "list_acronym": "TL",
            "display_name": "TestOrg",
            "category": "Viruses",
            "genus": "MyGenus",
            "species": "MySpecies",
        }
    ],
    "children": [
        {"child_taxid": "200", "controlled_taxid": "100", "child_name": "A Child"},
    ],
}


def _setup_regulated_state():
    """Populate module state directly (no filesystem) for unit-level query tests."""
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
        pd.DataFrame(
            [
                {
                    "list_acronym": "L1",
                    "list_item": "A.1: Flu A",
                    "tax_id": "11320",
                    "display_name": "Flu A",
                    "category": "Viruses",
                    "species": "FluSpecies",
                    "genus": "FluGenus",
                }
            ]
        )
    )
    ld.add_child_lut_data(
        pd.DataFrame(
            [
                {
                    "child_taxid": "99999",
                    "controlled_taxid": "11320",
                    "child_name": "Flu A Variant X",
                },
            ]
        )
    )
    tidy_control_list_data()


# ---------------------------------------------------------------------------
# import_data
# ---------------------------------------------------------------------------


def test_import_data_valid_path(tmp_path):
    _build_db(tmp_path, [_DEFAULT_LIST])
    assert import_data(str(tmp_path)) is True
    assert "TL" in ld.CONTROL_LISTS


def test_import_data_invalid_path():
    assert import_data("/nonexistent/path/to/nowhere") is False


def test_import_data_empty_directory(tmp_path):
    """A directory with region_definitions but no list folders yields False."""
    d = tmp_path / "empty"
    d.mkdir()
    (d / "region_definitions.json").write_text("[]")
    assert import_data(str(d)) is False


def test_import_data_nested_directories(tmp_path):
    """Lists inside nested subdirectories should be found via recursive search."""
    nested = tmp_path / "level1" / "level2"
    nested.mkdir(parents=True)
    (tmp_path / "region_definitions.json").write_text("[]")

    pd.DataFrame(
        [
            {
                "list_name": "TL",
                "list_acronym": "TL",
                "list_url": "url",
                "region_name": "NZ",
                "region_code": "NZ",
                "use": "EXPORT",
            }
        ]
    ).to_csv(str(nested / "list_info.csv"), index=False)
    pd.DataFrame(
        [
            {
                "tax_id": "100",
                "list_acronym": "TL",
                "display_name": "Org",
                "category": "Viruses",
                "species": "MySpecies",
                "genus": "MyGenus",
            }
        ]
    ).to_csv(str(nested / "controlled_taxids.csv"), index=False)
    pd.DataFrame(
        [
            {
                "child_taxid": "1",
                "controlled_taxid": "100",
                "child_name": "Org-like Borg Virus",
            }
        ]
    ).to_csv(str(nested / "children_of_controlled_taxids.csv"), index=False)

    assert import_data(str(tmp_path)) is True
    assert "TL" in ld.CONTROL_LISTS


def test_import_data_with_regional_context(tmp_path):
    _build_db(tmp_path, [_DEFAULT_LIST])
    assert import_data(str(tmp_path), regional_context=["NZ"]) is True
    assert ld.CONTROL_LISTS["TL"].status == ListMode.COMPLIANCE


def test_import_data_regional_context_non_matching(tmp_path):
    """Lists outside the regional context should be CONDITIONAL_NUM."""
    _build_db(tmp_path, [_DEFAULT_LIST])
    assert import_data(str(tmp_path), regional_context=["AU"]) is True
    assert ld.CONTROL_LISTS["TL"].status == ListMode.CONDITIONAL_NUM


# ---------------------------------------------------------------------------
# is_regulated
# ---------------------------------------------------------------------------


def test_is_regulated_true():
    _setup_regulated_state()
    assert is_regulated("11320") is True


def test_is_regulated_false():
    _setup_regulated_state()
    assert is_regulated("77777") is False


def test_is_regulated_via_child_mapping():
    """A child TaxID should be regulated through its parent mapping."""
    _setup_regulated_state()
    assert is_regulated(99999) is True, "TaxID passed as number not working."
    assert is_regulated("99999") is True, "TaxID passed as string not working."


# ---------------------------------------------------------------------------
# get_regulation
# ---------------------the------------------------------------------------------


def test_get_regulation_direct_hit():
    _setup_regulated_state()
    outputs = get_regulation("11320")
    assert len(outputs) == 1
    assert outputs[0].name == "Flu A"
    assert outputs[0].list == "L1"
    assert outputs[0].is_child is False


def test_get_regulation_child_hit():
    _setup_regulated_state()
    outputs = get_regulation("99999")
    assert len(outputs) == 1
    assert outputs[0].is_child is True


def test_get_regulation_no_hit():
    _setup_regulated_state()
    outputs = get_regulation("77777")
    assert len(outputs) == 0
    assert len(outputs) == 0


def test_get_regulation_multiple_lists():
    """A TaxID on two lists should produce two output entries."""
    ld.add_control_list(
        ControlList(
            "L1",
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
            "L2",
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
                {
                    "list_acronym": "L1",
                    "list_item": "A.1. Org Virus",
                    "tax_id": "100",
                    "display_name": "Org",
                    "category": "Viruses",
                    "species": "MySpecies",
                    "genus": "MyGenus",
                },
                {
                    "list_acronym": "L2",
                    "list_item": "C.5. Org Virus",
                    "tax_id": "100",
                    "display_name": "Org",
                    "category": "Viruses",
                    "species": "MySpecies",
                    "genus": "MyGenus",
                },
            ]
        )
    )

    tidy_control_list_data()

    outputs = get_regulation("100")
    assert len(outputs) == 2
    assert {o.list for o in outputs} == {"L1", "L2"}


# ---------------------------------------------------------------------------
# get_control_lists
# ---------------------------------------------------------------------------


def test_get_control_lists_all():
    ld.add_control_list(
        ControlList(
            "L1",
            "Lista Prima",
            "L1",
            "u1",
            Region("NZ", "NZ"),
            ListMode.COMPLIANCE,
            "EXPORT",
        )
    )
    ld.add_control_list(
        ControlList(
            "L2",
            "Lista Secunda",
            "L2",
            "u2",
            Region("AU", "AU"),
            ListMode.COMPLIANCE,
            "EXPORT",
        )
    )
    result = get_control_lists()
    assert isinstance(result, list)
    assert len(result) == 2


def test_get_control_lists_specific():
    cl = ControlList(
        "L1",
        "Lista Prima",
        "L1",
        "url",
        Region("NZ", "NZ"),
        ListMode.COMPLIANCE,
        "EXPORT",
    )
    ld.add_control_list(cl)
    result = get_control_lists("L1")
    assert result == cl


def test_get_control_lists_nonexistent():
    assert get_control_lists("NOPE") is None


def test_get_control_lists_empty():
    result = get_control_lists()
    assert result == []
