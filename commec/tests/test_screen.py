
import os
import json
from dataclasses import asdict
import textwrap
from unittest.mock import patch
import pytest
import pandas as pd
from commec.screen import run, ScreenArgumentParser, add_args
from commec.config.json_io import get_screen_data_from_json, encode_screen_data_to_json
from commec.utils.json_html_output import generate_html_from_screen_data
from commec.config.result import ScreenResult, ScreenStep, ScreenStatus

from commec.tests.screen_factory import ScreenTesterFactory

DATABASE_DIRECTORY = os.path.join(os.path.dirname(__file__), "test_dbs/")

def cheat_taxonomy_info(
    blast: pd.DataFrame,
    _regulated_taxids: list[str],
    _vaccine_taxids: list[str],
    _db_path: str | os.PathLike,
    _threads: int,
):
    """ 
    Used to override get_taxonomic_labels()
    Everything is a virus, and 
    anything with 12345 taxid is regulated.
    """
    blast["regulated"] = blast["subject tax ids"] == 12345
    blast["superkingdom"] = "Viruses"
    blast["phylum"] = "phylum"
    blast["genus"] = "genus"
    blast["species"] = "species"
    return blast

def print_tmp_path_contents(tmp_path):
    """ 
    Debug print, to see what was created within tmp_path
    """
    print(f"Contents of {tmp_path}:")
    for path in tmp_path.rglob("*"):  # Recursively list all files and directories
        print(path.relative_to(tmp_path), "->", "DIR" if path.is_dir() else "FILE")

def test_functional_screen(tmp_path, request):
    """
    Rewrite of the functional test using ScreenTestFactory
    """
    # File directories
    html_output_path = os.path.join(os.path.dirname(__file__), "test_data/functional")
    json_output_path = tmp_path / "functional.output.json"
    expected_json_output_path = os.path.join(os.path.dirname(__file__), "test_data/functional.json")

    # Screen Test Factory
    functional_test = ScreenTesterFactory("functional", tmp_path)
    functional_test.add_query("FCTEST1", 600)

    #Biorisk
    functional_test.add_hit(ScreenStep.BIORISK, "FCTEST1", 7, 95, "Toxin1", "ShouldntClear", 11111, description="LargeAreaFlag", score = 500, regulated = True, superkingdom = "Viruses", species = "horriblus")
    functional_test.add_hit(ScreenStep.BIORISK, "FCTEST1", 34, 65, "Toxin1", "ShouldntClear", 22222, description="SmallImportantFlag", score = 1000, regulated = True, superkingdom = "Viruses", species = "extra-horriblus")
    functional_test.add_hit(ScreenStep.BIORISK, "FCTEST1", 49, 80, "Toxin1", "ShouldTrim", 33333, description="SmallUnimportantTRIM", score = 100, regulated = True, superkingdom = "Viruses", species = "unimporticus")
    functional_test.add_hit(ScreenStep.BIORISK, "FCTEST1", 109, 191, "Toxin2", "ShouldWarn", 22222, description="WarningExample",score = 1000, regulated = False, superkingdom = "Viruses", species = "extra-horriblus-factor")
    functional_test.add_hit(ScreenStep.BIORISK, "FCTEST1", 593, 505, "Toxin3", "ShouldWarn", 11111, description="ReverseExample", score = 500, regulated = False, superkingdom = "Viruses", species = "horriblus-factor")
    # Protein Taxonomy
    functional_test.add_hit(ScreenStep.TAXONOMY_AA, "FCTEST1", 320, 380, "ShouldntClear", "NR_HIT_FLAG1", "12345", regulated = True, superkingdom = "Viruses", species = "regulaticus")
    functional_test.add_hit(ScreenStep.TAXONOMY_AA, "FCTEST1", 410, 490, "ShouldClearBySynBio", "NR_HIT_FLAG2", "12345", regulated = True, superkingdom = "Viruses", species = "regulaticus")
    functional_test.add_hit(ScreenStep.TAXONOMY_AA, "FCTEST1", 410, 500, "ShouldntClear", "NR_HIT_FLAG3", "12345", regulated = True, superkingdom = "Viruses", species = "regulaticus")
    functional_test.add_hit(ScreenStep.TAXONOMY_AA, "FCTEST1", 310, 370, "ShouldClear", "NR_HIT_FLAG4", "12346", regulated = True, superkingdom = "Viruses", species = "fine-icus")
    functional_test.add_hit(ScreenStep.TAXONOMY_AA, "FCTEST1", 340, 390, "ShouldMixedReg", "NR_HIT_MIXED1", "12347", regulated = True, superkingdom = "Viruses", species = "danger-poop")
    functional_test.add_hit(ScreenStep.TAXONOMY_AA, "FCTEST1", 340, 390, "ShouldMixednonReg", "NR_HIT_MIXED2", "12348", regulated = False, superkingdom = "Bacteria", species = "cute-happy-bacter")
    functional_test.add_hit(ScreenStep.TAXONOMY_AA, "FCTEST1", 340, 390, "ShouldMixedNonReg", "NR_HIT_MIXED3", "12349", regulated = False, superkingdom = "Bacteria", species = "poopicus")
    # Nucleotide Taxonomy
    functional_test.add_hit(ScreenStep.TAXONOMY_NT, "FCTEST1", 220, 280, "SUBJECT", "NT_HIT_FLAG1", "12345", regulated = True, superkingdom = "Viruses")
    functional_test.add_hit(ScreenStep.TAXONOMY_NT, "FCTEST1", 110, 190, "SUBJECT", "NT_HIT_FLAG2", "12345", regulated = True, superkingdom = "Viruses")
    functional_test.add_hit(ScreenStep.TAXONOMY_NT, "FCTEST1", 110, 200, "SUBJECT", "NT_HIT_FLAG3", "12345", regulated = True, superkingdom = "Viruses")
    functional_test.add_hit(ScreenStep.TAXONOMY_NT, "FCTEST1", 310, 390, "Main", "NT_HIT_MIXED", "12345", regulated = True, superkingdom = "Viruses")
    functional_test.add_hit(ScreenStep.TAXONOMY_NT, "FCTEST1", 310, 390, "NonRegMixedWithMain", "NT_HIT_MIXED2", "12346", regulated = False, superkingdom = "Bacteria")
    # Low Concern
    functional_test.add_hit(ScreenStep.LOW_CONCERN_PROTEIN, "FCTEST1", 202, 370, "Benign1", "Benign1", description = "BenignHMMClear")
    functional_test.add_hit(ScreenStep.LOW_CONCERN_RNA, "FCTEST1", 50, 150, "BENIGNRNA", "12346", description = "BenignCMTestOutput")
    functional_test.add_hit(ScreenStep.LOW_CONCERN_DNA, "FCTEST1", 410, 480, "BENIGNSYNBIO", "210", description = "BenignBlastClear")

    result = functional_test.run()

    # If we are writing exemplare data, do it in raw, to test the json_io simultaneously.
    gen_examples = request.config.getoption("--gen-examples")
    if gen_examples:
        encode_screen_data_to_json(result, expected_json_output_path)

    # Test results vs expected results.
    expected_screen_result : ScreenResult = get_screen_data_from_json(expected_json_output_path)
    actual_screen_result : ScreenResult = get_screen_data_from_json(json_output_path)
    sanitize_for_test(expected_screen_result)
    sanitize_for_test(actual_screen_result)

    # Generates .gitignored functional.html for quick human comparison.
    generate_html_from_screen_data(actual_screen_result, html_output_path)

    # Convert both original and retrieved data to dictionaries and compare
    assert asdict(expected_screen_result) == asdict(actual_screen_result), (
        f"Functional test does not match predicted output, fix code,"
        f" or if new output is expected, run with --gen-examples\n"
        f"Test JSON Reference data: \n{asdict(actual_screen_result)}\n"
        f"Test JSON output data: \n{asdict(expected_screen_result)}"
    )

def sanitize_for_test(screen_result: ScreenResult):
    """
    Remove arbitrary changes to the JSON that may arise during testing but are not
    relevant and should not be compared or versioned. 
    """
    # Runtime for pytest may be different
    screen_result.commec_info.time_taken = None
    screen_result.commec_info.date_run = None

    # All search tool versions etc may change.
    screen_result.commec_info.search_tool_info = None

    # Pytest increments the filename version, so ignore the input file.
    screen_result.query_info.file = "/test_placeholder/"

def test_screen_factory(tmp_path):
    my_factory = ScreenTesterFactory("test_01", tmp_path)
    my_factory.add_query("query_01", 500)
    my_factory.add_hit(ScreenStep.BIORISK, "query_01", 99, 399, "bad_risk", "BR500", 200, regulated = True)
    my_factory.add_hit(ScreenStep.TAXONOMY_AA, "query_01", 100, 300, "reg_gene", "ACC500", 500, "imaginary_species", regulated = True)
    my_factory.add_hit(ScreenStep.TAXONOMY_NT, "query_01", 1, 90, "reg_gene2", "ACC501", 501, "imaginary_species", regulated = True)
    my_factory.add_hit(ScreenStep.LOW_CONCERN_PROTEIN, "query_01", 97, 402, "safe_protein", "SF21YKN", 256, "safeicius")
    my_factory.add_hit(ScreenStep.LOW_CONCERN_DNA, "query_01", 97, 402, "safe_dna", "SF22YKN", 256, "safeicius")
    my_factory.add_hit(ScreenStep.LOW_CONCERN_RNA, "query_01", 97, 402, "safe_rna", "SF23YKN", 256, "safeicius")
    result = my_factory.run()
    assert result.queries["query_01"].status.screen_status == ScreenStatus.FLAG
    assert result.queries["query_01"].status.biorisk == ScreenStatus.FLAG
    assert result.queries["query_01"].status.protein_taxonomy == ScreenStatus.CLEARED_FLAG
    assert result.queries["query_01"].status.nucleotide_taxonomy == ScreenStatus.FLAG
    assert result.queries["query_01"].status.low_concern == ScreenStatus.FLAG

def test_different_regions(tmp_path):
    """
    Creates a single hit, with many regions, then a single clear on one of those regions.
    The low concern hit should only clear a single region, and not the whole query.
    Final result should be FLAG. 
    Tests that hits with multiple regions are being cleared correctly.
    """
    screen_test = ScreenTesterFactory("repeating_taxonomy", tmp_path)
    screen_test.add_query("query1",1000)
    screen_test.add_hit(ScreenStep.TAXONOMY_AA, "query1", 30, 90, "RegRepeat", "RR55", 500, regulated=True)
    screen_test.add_hit(ScreenStep.TAXONOMY_AA, "query1", 100, 170, "RegRepeat", "RR55", 500, regulated=True)
    screen_test.add_hit(ScreenStep.TAXONOMY_AA, "query1", 190, 260, "RegRepeat", "RR55", 500, regulated=True)
    screen_test.add_hit(ScreenStep.TAXONOMY_AA, "query1", 300, 390, "RegRepeat", "RR55", 500, regulated=True)
    screen_test.add_hit(ScreenStep.TAXONOMY_AA, "query1", 400, 750, "RegRepeat", "RR55", 500, regulated=True)
    screen_test.add_hit(ScreenStep.LOW_CONCERN_PROTEIN, "query1", 400, 750, "ClearProtein", "RR55CLEAR", 500)
    result = screen_test.run()

    num_hits = len(result.queries["query1"].hits)
    num_regions = len(result.queries["query1"].hits["RR55"].ranges)
    status = result.queries["query1"].status.screen_status
    assert  num_hits == 2, f"Expected two hits, got {num_hits}."
    assert  num_regions == 5, (f"Number of ranges [{num_regions}] in hit `RR55` for "
                            " `query1` not equal to expected number (5).")
    assert status == ScreenStatus.FLAG, ("Expected status is to FLAG, current status"
                                        f" is {status}, likely multiple region clearing issue.")
