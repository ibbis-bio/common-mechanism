#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Functional tests for commec screen.
"""

import os
import json
from dataclasses import asdict
from commec.config.json_io import get_screen_data_from_json, encode_screen_data_to_json
from commec.utils.json_html_output import generate_html_from_screen_data
from commec.config.result import ScreenResult, ScreenStep, ScreenStatus
from commec.tests.screen_factory import ScreenTesterFactory

DATABASE_DIRECTORY = os.path.join(os.path.dirname(__file__), "test_dbs/")


def sanitize_for_test(screen_result: ScreenResult):
    """
    Remove arbitrary changes to the JSON that may arise during testing but are not
    relevant and should not be compared or versioned.
    """
    # Runtime for pytest may be different
    screen_result.commec_info.time_taken = None
    screen_result.commec_info.date_run = None

    # All search tool versions etc may change.
    screen_result.database_info = None

    # Pytest increments the filename version, so ignore the input file.
    screen_result.query_info.file = "/test_placeholder/"


def test_functional_screen(tmp_path, request):
    """
    Rewrite of the functional test using ScreenTestFactory
    """
    # File directories
    html_output_path = os.path.join(os.path.dirname(__file__), "test_data/functional")
    json_output_path = tmp_path / "functional.output.json"
    expected_json_output_path = os.path.join(
        os.path.dirname(__file__), "test_data/functional.json"
    )

    # Screen Test Factory
    functional_test = ScreenTesterFactory("functional", tmp_path)
    functional_test.add_query("FCTEST1", 600)

    # Biorisk
    functional_test.add_hit(
        ScreenStep.BIORISK,
        "FCTEST1",
        7,
        95,
        "Toxin1a",
        "ShouldntClear",
        11111,
        description="ShouldntClear LargeAreaFlag",
        score=500,
        regulated=True,
        superkingdom="Viruses",
        species="horriblus",
        genus="horribluses",
        evalue=1e-21,
    )
    functional_test.add_hit(
        ScreenStep.BIORISK,
        "FCTEST1",
        34,
        65,
        "Toxin1b",
        "ShouldntClearSmall",
        22222,
        description="ShouldntClear SmallImportantFlag",
        score=1000,
        regulated=True,
        superkingdom="Viruses",
        species="extra-horriblus",
        genus="horribluses",
        evalue=1e-22,
    )
    functional_test.add_hit(
        ScreenStep.BIORISK,
        "FCTEST1",
        49,
        80,
        "Toxin1c",
        "ShouldTrim",
        33333,
        description="ShouldTrim SmallUnimportantTRIM",
        score=100,
        regulated=True,
        superkingdom="Viruses",
        species="unimporticus",
        genus="pseudohorribluses",
        evalue=1e-23,
    )
    functional_test.add_hit(
        ScreenStep.BIORISK,
        "FCTEST1",
        109,
        191,
        "Toxin2",
        "ShouldWarn",
        22222,
        description="WarningExample",
        score=1000,
        regulated=False,
        superkingdom="Viruses",
        species="extra-horriblus-factor",
        genus="horribluses",
        evalue=1e-24,
    )
    functional_test.add_hit(
        ScreenStep.BIORISK,
        "FCTEST1",
        593,
        505,
        "Toxin3",
        "ShouldWarn2",
        11111,
        description="ReverseExample",
        score=500,
        regulated=False,
        superkingdom="Viruses",
        species="horriblus-factor",
        genus="horribluses",
        evalue=1e-25,
    )
    # Protein Taxonomy
    functional_test.add_hit(
        ScreenStep.TAXONOMY_AA,
        "FCTEST1",
        320,
        380,
        "NR_HIT_ShouldntClear",
        "NR_HIT_FLAG1",
        "12345",
        regulated=True,
        superkingdom="Viruses",
        species="regulaticus",
        genus="orthoregulatidae",
        evalue=0.06,
    )
    functional_test.add_hit(
        ScreenStep.TAXONOMY_AA,
        "FCTEST1",
        410,
        490,
        "ShouldClearBySynBio",
        "NR_HIT_FLAG2",
        "12347",
        regulated=True,
        superkingdom="Viruses",
        species="regulaticus",
        genus="orthoregulatidae",
        evalue=0.07,
    )
    functional_test.add_hit(
        ScreenStep.TAXONOMY_AA,
        "FCTEST1",
        410,
        500,
        "NR_HIT2_ShouldntClear",
        "NR_HIT_FLAG3",
        "12345",
        regulated=True,
        superkingdom="Viruses",
        species="regulaticus",
        genus="orthoregulatidae",
        evalue=0.08,
    )
    functional_test.add_hit(
        ScreenStep.TAXONOMY_AA,
        "FCTEST1",
        310,
        370,
        "ShouldClear",
        "NR_HIT_FLAG4",
        "12346",
        regulated=True,
        superkingdom="Viruses",
        species="fine-icus",
        genus="orthoderegulatidae",
        evalue=0.09,
    )
    functional_test.add_hit(
        ScreenStep.TAXONOMY_AA,
        "FCTEST1",
        340,
        390,
        "ShouldMixedReg1",
        "NR_HIT_MIXED1",
        "12347",
        regulated=True,
        superkingdom="Viruses",
        species="danger-poop",
        genus="faecia",
        evalue=0.10,
    )
    functional_test.add_hit(
        ScreenStep.TAXONOMY_AA,
        "FCTEST1",
        340,
        390,
        "ShouldMixednonReg2",
        "NR_HIT_MIXED2",
        "12348",
        regulated=False,
        superkingdom="Bacteria",
        species="cute-happy-bacter",
        genus="happicillus",
        evalue=0.11,
    )
    functional_test.add_hit(
        ScreenStep.TAXONOMY_AA,
        "FCTEST1",
        340,
        390,
        "ShouldMixedNonReg3",
        "NR_HIT_MIXED3",
        "12349",
        regulated=False,
        superkingdom="Bacteria",
        species="poopicus",
        genus="faecia",
        evalue=0.12,
    )
    # Nucleotide Taxonomy
    functional_test.add_hit(
        ScreenStep.TAXONOMY_NT,
        "FCTEST1",
        220,
        280,
        "SUBJECT_1",
        "NT_HIT_FLAG1",
        "22345",
        description="Should clear from RNA (its on the cusp)",
        regulated=True,
        superkingdom="Viruses",
        species="Viria Novus",
        genus="orthonovaregulatidae",
        evalue=0.13,
    )
    functional_test.add_hit(
        ScreenStep.TAXONOMY_NT,
        "FCTEST1",
        110,
        190,
        "SUBJECT_2",
        "NT_HIT_FLAG2",
        "22345",
        regulated=True,
        superkingdom="Viruses",
        species="Viria Novus",
        genus="orthonovaregulatidae",
        evalue=0.14,
    )
    functional_test.add_hit(
        ScreenStep.TAXONOMY_NT,
        "FCTEST1",
        110,
        200,
        "SUBJECT_3",
        "NT_HIT_FLAG3",
        "12350",
        description="Should Flag",
        regulated=True,
        superkingdom="Viruses",
        species="Viria Novia blue beads",
        genus="orthonovaregulatidae",
        evalue=0.15,
    )
    functional_test.add_hit(
        ScreenStep.TAXONOMY_NT,
        "FCTEST1",
        310,
        390,
        "Main",
        "NT_HIT_MIXED",
        "12345",
        regulated=True,
        superkingdom="Viruses",
        species="Viria Novus",
        genus="orthonovaregulatidae",
        evalue=0.16,
    )
    functional_test.add_hit(
        ScreenStep.TAXONOMY_NT,
        "FCTEST1",
        310,
        390,
        "NonRegMixedWithMain",
        "NT_HIT_MIXED2",
        "12348",
        regulated=False,
        superkingdom="Bacteria",
        species="SafetyGreenius Virus",
        genus="Greenocollii",
        evalue=0.17,
    )
    # Low Concern
    functional_test.add_hit(
        ScreenStep.LOW_CONCERN_PROTEIN,
        "FCTEST1",
        202,
        370,
        "BENIGNPROT",
        "Benign1",
        description="BenignHMMClear",
        evalue=1e-26,
    )
    functional_test.add_hit(
        ScreenStep.LOW_CONCERN_RNA,
        "FCTEST1",
        50,
        150,
        "BENIGNRNA",
        "211",
        description="BenignCMTestOutput",
        evalue=1e-27,
    )
    functional_test.add_hit(
        ScreenStep.LOW_CONCERN_DNA,
        "FCTEST1",
        410,
        480,
        "BENIGNSYNBIO",
        "210",
        description="BenignBlastClear",
        evalue=1e-28,
    )

    result = functional_test.run()

    # If we are writing exemplare data, do it in raw, to test the json_io simultaneously.
    gen_examples = request.config.getoption("--gen-examples")
    if gen_examples:
        encode_screen_data_to_json(result, expected_json_output_path)

    # Test results vs expected results.
    expected_screen_result: ScreenResult = get_screen_data_from_json(
        expected_json_output_path
    )
    actual_screen_result: ScreenResult = get_screen_data_from_json(json_output_path)
    sanitize_for_test(expected_screen_result)
    sanitize_for_test(actual_screen_result)

    # Generates .gitignored functional.html for quick human comparison.
    generate_html_from_screen_data(actual_screen_result, html_output_path)

    # Convert both original and retrieved data to dictionaries and compare
    assert asdict(expected_screen_result) == asdict(actual_screen_result), (
        "Functional test does not match predicted output, fix code,"
        " or if new output is expected, run with --gen-examples\n"
    )


def test_screen_factory(tmp_path):
    my_factory = ScreenTesterFactory("test_01", tmp_path)
    my_factory.add_query("query01", 500)
    my_factory.add_hit(
        ScreenStep.BIORISK, "query01", 99, 399, "bad_risk", "BR500", 200, regulated=True
    )
    my_factory.add_hit(
        ScreenStep.TAXONOMY_AA,
        "query01",
        100,
        300,
        "reg_gene",
        "ACC500",
        500,
        "imaginary_species",
        regulated=True,
    )
    my_factory.add_hit(
        ScreenStep.TAXONOMY_NT,
        "query01",
        1,
        90,
        "reg_gene2",
        "ACC501",
        501,
        "imaginary_species",
        regulated=True,
    )
    my_factory.add_hit(
        ScreenStep.LOW_CONCERN_PROTEIN,
        "query01",
        97,
        402,
        "safe_protein",
        "SF21YKN",
        256,
        "safeicius",
    )
    my_factory.add_hit(
        ScreenStep.LOW_CONCERN_DNA,
        "query01",
        97,
        402,
        "safe_dna",
        "SF22YKN",
        256,
        "safeicius",
    )
    my_factory.add_hit(
        ScreenStep.LOW_CONCERN_RNA,
        "query01",
        97,
        402,
        "safe_rna",
        "SF23YKN",
        256,
        "safeicius",
    )
    result = my_factory.run()
    assert result.queries["query01"].status.screen_status == ScreenStatus.FLAG
    assert result.queries["query01"].status.biorisk == ScreenStatus.FLAG
    assert (
        result.queries["query01"].status.protein_taxonomy == ScreenStatus.CLEARED_FLAG
    )
    assert result.queries["query01"].status.nucleotide_taxonomy == ScreenStatus.FLAG
    assert result.queries["query01"].status.low_concern == ScreenStatus.FLAG


def test_low_concern_multiquery(tmp_path):
    """
    Checks that the low concern hits are correctly accessed when multiple queries are used, and
    the low concern databases have multiple hits across all queries.
    """
    screen_test = ScreenTesterFactory("lowconcern_multiquerycheck", tmp_path)
    screen_test.add_query("query1", 1000)
    screen_test.add_query("query2", 1000)
    screen_test.add_query("query3", 1000)
    screen_test.add_hit(
        ScreenStep.TAXONOMY_AA,
        "query1",
        400,
        600,
        "ClearMe",
        "RR55",
        500,
        regulated=True,
    )
    screen_test.add_hit(
        ScreenStep.TAXONOMY_AA,
        "query2",
        400,
        600,
        "ClearMe",
        "RR55",
        500,
        regulated=True,
    )
    screen_test.add_hit(
        ScreenStep.TAXONOMY_AA,
        "query3",
        400,
        600,
        "ClearMe",
        "RR55",
        500,
        regulated=True,
    )
    screen_test.add_hit(
        ScreenStep.LOW_CONCERN_RNA, "query1", 100, 900, "ClearProtein", "RR55CLEAR", 500
    )
    screen_test.add_hit(
        ScreenStep.LOW_CONCERN_RNA, "query2", 100, 900, "ClearProtein", "RR55CLEAR", 500
    )
    screen_test.add_hit(
        ScreenStep.LOW_CONCERN_RNA, "query3", 100, 900, "ClearProtein", "RR55CLEAR", 500
    )
    result = screen_test.run()

    assert result.queries["query1"].status.screen_status == ScreenStatus.CLEARED_FLAG
    assert result.queries["query2"].status.screen_status == ScreenStatus.CLEARED_FLAG
    assert result.queries["query3"].status.screen_status == ScreenStatus.CLEARED_FLAG


def test_different_regions(tmp_path):
    """
    Creates a single hit, with many regions, then a single clear on one of those regions.
    The low concern hit should only clear a single region, and not the whole query.
    Final result should be FLAG.
    Tests that hits with multiple regions are being cleared correctly.
    """
    screen_test = ScreenTesterFactory("repeating_taxonomy", tmp_path)
    screen_test.add_query("query1", 1000)
    screen_test.add_hit(
        ScreenStep.TAXONOMY_AA,
        "query1",
        30,
        90,
        "RegRepeat",
        "12345",
        500,
        regulated=True,
    )
    screen_test.add_hit(
        ScreenStep.TAXONOMY_AA,
        "query1",
        100,
        189,
        "RegRepeat",
        "12345",
        500,
        regulated=True,
    )
    screen_test.add_hit(
        ScreenStep.TAXONOMY_AA,
        "query1",
        190,
        260,
        "RegRepeat",
        "12345",
        500,
        regulated=True,
    )
    screen_test.add_hit(
        ScreenStep.TAXONOMY_AA,
        "query1",
        200,
        380,
        "RegRepeat",
        "12345",
        500,
        regulated=True,
    )
    screen_test.add_hit(
        ScreenStep.TAXONOMY_AA,
        "query1",
        400,
        750,
        "RegRepeat",
        "12345",
        500,
        regulated=True,
    )
    screen_test.add_hit(
        ScreenStep.LOW_CONCERN_PROTEIN, "query1", 400, 750, "ClearProtein", "12346", 500
    )
    result = screen_test.run()

    # encode_screen_data_to_json(result, "../test_output.json")

    num_hits = len(result.queries["query1"].hits)
    status = result.queries["query1"].status.screen_status
    assert num_hits == 4, (
        f"Got {num_hits}, instead of 4 hits. Expected 3 Clusters for tax_id 500, and 1 benign hit...\n{json.dumps(asdict(result.queries['query1']), indent=2)}"
    )
    assert status == ScreenStatus.FLAG, (
        "Expected status is to FLAG, current status"
        f" is {status}, likely multiple region clearing issue."
    )


def test_no_hits_outcomes(tmp_path):
    """
    Simply tests what occurs when commec is run and finds no hits, with and without --skip-tx
    """
    screen_test = ScreenTesterFactory("no_hits_test", tmp_path)
    screen_test.add_query("query_not_hits", 1000)
    result = screen_test.run()
    assert result.queries["query_not_hits"].status.screen_status == ScreenStatus.PASS, (
        f"Query that had no hits, didn't have correct screen status.: {result.queries['query_not_hits'].status.screen_status}"
    )

    screen_test_skip_tx = ScreenTesterFactory("no_hits_test_skip_tx", tmp_path)
    screen_test_skip_tx.add_query("query_not_hits", 1000)
    result = screen_test_skip_tx.run("--skip-tx")
    assert (
        result.queries["query_not_hits"].status.screen_status
        == ScreenStatus.PASS_SKIP_TX
    ), (
        f"Query that skipped Tx, but had no hits, didn't have correct screen status.: {result.queries['query_not_hits'].status.screen_status}"
    )
