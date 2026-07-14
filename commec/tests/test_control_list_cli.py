"""Tests for commec.control_list.cli — CLI argument handling and output formatting."""

import argparse
import json
import os

import pytest
import pandas as pd

import commec.control_list.list_data as ld
import commec.control_list.region as region_module
from commec.control_list.containers import (
    ControlList,
    ControlListContext,
    ControlListOutput,
    ListMode,
    Region,
)
from commec.control_list.cli import (
    add_args,
    format_control_list_annotation,
    format_control_lists,
    generate_output_summary_csv,
)
from commec.control_list.control_list import run
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


def _setup_lists_with_annotations():
    """Populate module state with two lists and a few annotations."""
    ld.add_control_list(ControlList(
        "Export Controls", "Moderationes Exportationis", "EC", "http://ec.gov",
        Region("New Zealand", "NZ"),
        ListMode.COMPLIANCE, "EXPORT"))
    ld.add_control_list(ControlList(
        "Pathogen List", "Charta Pathogenorum", "PL", "http://pl.gov",
        Region("Australia", "AU"),
        ListMode.CONDITIONAL_NUM, "PATHGN"))
    ld.add_control_list_annotations(pd.DataFrame([
        {"list_acronym": "EC", "list_item": "A.1. Virus A", "tax_id": "100", "display_name": "Virus A",
         "category": "Viruses", "species" : "Aicus Virus", "genus" : "BadVirusGenus"},
        {"list_acronym": "EC", "list_item": "A.2. Virus B", "tax_id": "200", "display_name": "Virus B",
         "category": "Viruses", "species" : "Bicus Virus", "genus" : "BadVirusGenus"},
        {"list_acronym": "PL", "list_item": "G.10. Bacterium X", "tax_id": "300", "display_name": "Bacterium X",
         "category": "Bacteria", "species" : "Baccilus Xanthomouse", "genus" : "BadBacteriaGenus"},
    ]))
    tidy_control_list_data()


def _build_minimal_db(base_path):
    """Build the smallest valid control-list database directory."""
    os.makedirs(str(base_path), exist_ok=True)
    with open(os.path.join(str(base_path), "region_definitions.json"), "w") as f:
        json.dump([], f)
    list_dir = os.path.join(str(base_path), "testlist")
    os.makedirs(list_dir, exist_ok=True)
    pd.DataFrame([{
        "list_name": "TL", "list_acronym": "TL", "list_url": "url",
        "region_name": "NZ", "region_code": "NZ", "use": "EXPORT",
    }]).to_csv(os.path.join(list_dir, "list_info.csv"), index=False)
    pd.DataFrame([{
        "tax_id": "100", "list_acronym": "TL", "display_name": "Org",
        "category": "Viruses", "species" : "Orgonovirus Organa", "genus" : "Orthoorgonovirae",
    }]).to_csv(os.path.join(list_dir, "controlled_taxids.csv"), index=False)
    pd.DataFrame([{
        "child_taxid": "1", "controlled_taxid": "100", "child_name" : "Orgonovirus Sisterina"
    }]).to_csv(os.path.join(list_dir, "children_of_controlled_taxids.csv"), index=False)


# ---------------------------------------------------------------------------
# format_control_lists
# ---------------------------------------------------------------------------

def test_format_control_lists_table():
    _setup_lists_with_annotations()
    output = format_control_lists(verbosity=False)
    assert "EC" in output
    assert "PL" in output
    assert "Control List" in output
    assert "Acronym" in output


def test_format_control_lists_empty():
    """With no lists loaded, should return a string without raising."""
    output = format_control_lists(verbosity=False)
    assert isinstance(output, str)


def test_format_control_lists_verbose_empty():
    output = format_control_lists(verbosity=True)
    assert "The following Control Lists apply" in output
    assert "Total number of Taxid Relationships:0" in output


# ---------------------------------------------------------------------------
# format_control_list_annotation
# ---------------------------------------------------------------------------

def test_format_control_list_annotation_single():
    data = [ControlListOutput("Flu A", "Viruses", "EC", "flu_species", "flu_genus")]
    output = format_control_list_annotation(data)
    assert "viruses" in output # Its lowered in the string output.
    assert "Flu A" in output
    assert "EC" in output
    # Single entry should NOT have the plural header
    assert "Controlled by the following lists" not in output


def test_format_control_list_annotation_multiple():
    data = [
        ControlListOutput("Flu A", "Viruses", "EC"),
        ControlListOutput("Flu A", "Viruses", "PL"),
    ]
    output = format_control_list_annotation(data)
    assert "Controlled by the following lists" in output
    assert "EC" in output
    assert "PL" in output


def test_format_control_list_annotation_with_derived_from():
    data = [ControlListOutput("Flu A", "Viruses", "EC","","","Parent Organism",False,"Child Name")]
    output = format_control_list_annotation(data)
    assert "Parent Organism" in output
    assert "Child Name" not in output


# ---------------------------------------------------------------------------
# generate_output_summary_csv
# ---------------------------------------------------------------------------

def test_generate_output_summary_csv(tmp_path):
    _setup_lists_with_annotations()
    output_path = tmp_path / "summary.csv"
    generate_output_summary_csv(str(output_path))
    assert output_path.exists()
    df = pd.read_csv(output_path)
    # One-hot encoded list acronym columns should be present
    assert "EC" in df.columns
    assert "PL" in df.columns


def test_generate_output_summary_csv_suffix_correction(tmp_path):
    """A path without .csv suffix should be auto-corrected."""
    _setup_lists_with_annotations()
    generate_output_summary_csv(str(tmp_path / "summary"))
    assert (tmp_path / "summary.csv").exists()


def test_generate_output_summary_csv_row_count(tmp_path):
    """Each unique accession should produce one row in the output CSV."""
    _setup_lists_with_annotations()
    output_path = tmp_path / "out.csv"
    generate_output_summary_csv(str(output_path))
    df = pd.read_csv(output_path)
    # We set up 3 distinct taxids (100, 200, 300)
    assert len(df) == 3


# ---------------------------------------------------------------------------
# add_args
# ---------------------------------------------------------------------------

def test_add_args_returns_parser():
    parser = argparse.ArgumentParser()
    result = add_args(parser)
    assert result is parser


def test_add_args_list_flag():
    parser = argparse.ArgumentParser()
    add_args(parser)
    args = parser.parse_args(["--list"])
    assert args.showlists is True
    assert args.verbose is False


def test_add_args_accessions():
    parser = argparse.ArgumentParser()
    add_args(parser)
    args = parser.parse_args(["--accessions", "12345", "67890"])
    assert args.showtaxids == ["12345", "67890"]


def test_add_args_regions():
    parser = argparse.ArgumentParser()
    add_args(parser)
    args = parser.parse_args(["--regions", "NZ", "AU"])
    assert args.regions == ["NZ", "AU"]


def test_add_args_verbose():
    parser = argparse.ArgumentParser()
    add_args(parser)
    args = parser.parse_args(["--verbose", "--list"])
    assert args.verbose is True


def test_add_args_output_prefix():
    parser = argparse.ArgumentParser()
    add_args(parser)
    args = parser.parse_args(["--output_prefix", "myprefix", "--list"])
    assert args.output_prefix == "myprefix"


def test_add_args_databases_and_config_mutually_exclusive(tmp_path):
    """Providing both -d and -y should cause argparse to error."""
    parser = argparse.ArgumentParser()
    add_args(parser)
    with pytest.raises(SystemExit):
        parser.parse_args(["-d", str(tmp_path), "-y", str(tmp_path / "c.yaml")])


# ---------------------------------------------------------------------------
# run() — CLI entry-point
# ---------------------------------------------------------------------------

def test_run_cli_no_action():
    """run() returns 1 when neither --list nor --accessions or --output_prefix is specified."""
    args = argparse.Namespace(
        verbose=False, database_dir=None, yaml_file=None,
        regions=[], showlists=False, showtaxids=None, output_prefix="",
    )
    assert run(args) == 1


def test_run_cli_no_database():
    """run() returns 1 when an action is requested but no database is given."""
    args = argparse.Namespace(
        verbose=False, database_dir=None, yaml_file=None,
        regions=[], showlists=True, showtaxids=None, output_prefix="",
    )
    assert run(args) == 2


def test_run_cli_with_list_flag(tmp_path):
    db_path = tmp_path / "db"
    _build_minimal_db(db_path)
    args = argparse.Namespace(
        verbose=False, database_dir=str(db_path), yaml_file=None,
        regions=[], showlists=True, showtaxids=None, output_prefix="",
    )
    result = run(args)
    # Successful run returns None (no explicit return value)
    assert result is None


def test_run_cli_with_accessions(tmp_path):
    db_path = tmp_path / "db"
    _build_minimal_db(db_path)
    args = argparse.Namespace(
        verbose=False, database_dir=str(db_path), yaml_file=None,
        regions=[], showlists=False, showtaxids=["100"], output_prefix="",
    )
    result = run(args)
    assert result is None


def test_run_cli_with_output_prefix(tmp_path):
    db_path = tmp_path / "db"
    _build_minimal_db(db_path)
    output_prefix = str(tmp_path / "output_summary")
    args = argparse.Namespace(
        verbose=False, database_dir=str(db_path), yaml_file=None,
        regions=[], showlists=True, showtaxids=None, output_prefix=output_prefix,
    )
    run(args)
    assert (tmp_path / "output_summary.csv").exists()
