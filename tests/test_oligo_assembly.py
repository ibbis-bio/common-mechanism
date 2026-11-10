"""Verify source attribution without using private screening sequences."""

import copy

import pytest

from commec.oligo_assembly import assemble, combine, load_records, project_region


def fixture_data(start=111, end=150, hit_status="Flag", contig_status="Flag"):
    records = [
        {
            "id": f"q{i}",
            "original_id": "same_header",
            "header": "same_header",
            "record_index": i,
            "length": 100,
        }
        for i in range(3)
    ]
    direct = {
        "queries": {
            r["id"]: {"query": r["id"], "hits": [], "status": {"screen_status": "Pass"}}
            for r in records
        }
    }
    hit = {
        "name": "Synthetic test hit",
        "recommendation": {
            "status": hit_status,
            "from_step": "Protein Taxonomy Search",
        },
        "region": {"query_start": start, "query_end": end},
        "annotations": {},
    }
    assemblies = {
        "queries": {"a0": {"status": {"screen_status": contig_status}, "hits": [hit]}}
    }
    manifest = {
        "assembler": "oligraph",
        "version": "0.2.2",
        "min_overlap": 15,
        "matching": "exact",
        "method": "all",
        "records": records,
        "exclusions": [],
        "contigs": [
            {
                "id": "a0",
                "length": 240,
                "placements": [
                    {"id": "q0", "strand": "+", "start": 0, "end": 100},
                    {"id": "q1", "strand": "+", "start": 80, "end": 180},
                    {"id": "q2", "strand": "+", "start": 160, "end": 240},
                ],
            }
        ],
    }
    return direct, assemblies, manifest


def test_only_the_oligo_overlapping_the_concern_is_flagged():
    direct, assemblies, manifest = fixture_data()
    original = copy.deepcopy(direct)
    effective, summary = combine(direct, assemblies, manifest)
    assert [r["final_status"] for r in summary["oligos"]] == ["Pass", "Flag", "Pass"]
    assert summary["flagged_oligos"] == 1
    assert summary["pool_verdict"] is None
    assert direct == original
    assert (
        effective["queries"]["q1"]["assembly"]["evidence"][0]["hit"]["region"][
            "query_start"
        ]
        == 31
    )
    assert (
        effective["queries"]["q1"]["assembly"]["evidence"][0]["hit"]["region"][
            "query_end"
        ]
        == 70
    )


def test_shared_overlap_attributes_both_actual_sources_once():
    direct, assemblies, manifest = fixture_data(start=91, end=95)
    _, summary = combine(direct, assemblies, manifest)
    assert [r["final_status"] for r in summary["oligos"]] == ["Flag", "Flag", "Pass"]
    assert summary["flagged_oligos"] == 2


def test_direct_flag_is_preserved_when_assembly_has_no_active_concern():
    direct, assemblies, manifest = fixture_data(
        hit_status="Flag (Cleared)", contig_status="Pass"
    )
    direct["queries"]["q0"]["status"]["screen_status"] = "Flag"
    _, summary = combine(direct, assemblies, manifest)
    assert [r["final_status"] for r in summary["oligos"]] == ["Flag", "Pass", "Pass"]


def test_warning_and_short_query_status_are_not_silently_cleared():
    direct, assemblies, manifest = fixture_data(
        hit_status="Warning", contig_status="Warning"
    )
    direct["queries"]["q2"]["status"]["screen_status"] = "Skip (too short)"
    _, summary = combine(direct, assemblies, manifest)
    assert [r["final_status"] for r in summary["oligos"]] == [
        "Pass",
        "Warning",
        "Skip (too short)",
    ]
    assert summary["oligos"][2]["assembly_coverage"] == "screened"


def test_reverse_strand_and_half_open_boundaries():
    placement = {"start": 80, "end": 180, "strand": "-"}
    assert project_region(placement, 90, 110) == (70, 90)
    assert project_region(placement, 180, 200) is None
    assert project_region(placement, 60, 80) is None


def test_invalid_hit_coordinates_are_reported_without_flagging_every_source():
    direct, assemblies, manifest = fixture_data(start=0, end=900)
    _, summary = combine(direct, assemblies, manifest)
    assert summary["flagged_oligos"] == 0
    assert len(summary["unresolved_attributions"]) == 1
    assert all(row["final_status"] == "Incomplete" for row in summary["oligos"])


def test_lost_assembly_result_fails_the_accounting_check():
    direct, assemblies, manifest = fixture_data()
    assemblies["queries"].clear()
    with pytest.raises(ValueError, match="every contig"):
        combine(direct, assemblies, manifest)


def test_duplicate_headers_keep_distinct_record_identities(tmp_path):
    path = tmp_path / "input.fasta"
    path.write_text(">duplicate first\nACGT\n>duplicate second\nTGCA\n")
    records = load_records(path)
    assert records[0]["original_id"] == records[1]["original_id"]
    assert records[0]["id"] != records[1]["id"]
    assert records[0]["header"] != records[1]["header"]


def test_real_assembler_replays_a_benign_reverse_complement_path():
    pytest.importorskip("oligraph._oligraph")
    from Bio.Seq import Seq

    reference = "ACTGCTTACGAGCTGATCGTACCGATGCTAGCTTGACGTAGCCTATGACCTGAGCATCGATGAC"
    records = [
        {
            "id": "q0",
            "header": "left",
            "original_id": "left",
            "record_index": 0,
            "sequence": reference[:45],
        },
        {
            "id": "q1",
            "header": "right",
            "original_id": "right",
            "record_index": 1,
            "sequence": str(Seq(reference[25:]).reverse_complement()),
        },
    ]
    manifest, contigs = assemble(records)
    assert len(contigs) == 1
    assert contigs[0]["sequence"] in {
        reference,
        str(Seq(reference).reverse_complement()),
    }
    assert len(manifest["contigs"][0]["placements"]) == 2


def test_contig_error_keeps_known_flag_and_marks_other_sources_incomplete():
    direct, assemblies, manifest = fixture_data(contig_status="Error")
    _, summary = combine(direct, assemblies, manifest)
    assert [r["final_status"] for r in summary["oligos"]] == [
        "Incomplete",
        "Flag",
        "Incomplete",
    ]
    assert all(row["assembly_coverage"] == "incomplete" for row in summary["oligos"])
    assert summary["oligos"][1]["assembly_status"] == "Flag"
    assert summary["oligos"][0]["assembly_status"] == "No attributed hit"
    assert summary["flagged_assemblies"] == 1
    assert summary["assembly_statuses"] == {"a0": "Error"}


def test_direct_error_keeps_known_flag_evidence():
    direct, assemblies, manifest = fixture_data(
        hit_status="Flag (Cleared)", contig_status="Pass"
    )
    direct["queries"]["q0"]["status"]["screen_status"] = "Error"
    direct["queries"]["q0"]["hits"] = [{"recommendation": {"status": "Flag"}}]
    _, summary = combine(direct, assemblies, manifest)
    assert summary["oligos"][0]["final_status"] == "Error"
    assert summary["oligos"][0]["incomplete_screening"]
    assert summary["newly_flagged_oligos"] == 0


def test_known_flag_is_not_downgraded_by_incomplete_assembly():
    direct, assemblies, manifest = fixture_data(contig_status="Error")
    direct["queries"]["q0"]["status"]["screen_status"] = "Flag"
    _, summary = combine(direct, assemblies, manifest)
    assert summary["oligos"][0]["final_status"] == "Flag"
    assert summary["oligos"][1]["final_status"] == "Flag"
    assert all(row["assembly_coverage"] == "incomplete" for row in summary["oligos"])


def test_screened_contig_does_not_erase_another_incomplete_contig():
    direct, assemblies, manifest = fixture_data()
    second = copy.deepcopy(manifest["contigs"][0])
    second["id"] = "a1"
    manifest["contigs"].append(second)
    assemblies["queries"]["a0"]["status"]["screen_status"] = "Error"
    assemblies["queries"]["a1"] = copy.deepcopy(assemblies["queries"]["a0"])
    assemblies["queries"]["a1"]["status"]["screen_status"] = "Pass"
    _, summary = combine(direct, assemblies, manifest)
    assert all(row["assembly_coverage"] == "incomplete" for row in summary["oligos"])


def test_all_excluded_records_keep_reasons_without_assembly_failure():
    records = [{"id": "q0", "sequence": "N" * 30}]
    manifest, contigs = assemble(records)
    assert contigs == []
    assert manifest["exclusions"][0]["id"] == "q0"
    assert "error" not in manifest


def test_html_shows_duplicate_record_numbers_and_assembly_skip(tmp_path):
    from dataclasses import asdict

    from commec.config.result import ScreenResult
    from commec.oligo_assembly import render, write_json

    direct, assemblies, manifest = fixture_data(
        hit_status="Flag (Cleared)", contig_status="Skip (too short)"
    )
    effective, summary = combine(direct, assemblies, manifest)
    effective_path, output = tmp_path / "effective.json", tmp_path / "report.html"
    document = asdict(ScreenResult())
    document["queries"] = effective["queries"]
    document["assembly_info"] = effective["assembly_info"]
    write_json(effective_path, document)
    render(effective_path, summary, manifest, output)
    page = output.read_text()
    assert "same_header (record 1)" in page
    assert "same_header (record 2)" in page
    assert '"enabled": true' in page
    assert '"coverage": "incomplete"' in page
    assert "Assembly screening proof of concept" not in page


@pytest.mark.parametrize("contig_status", ["Skip (too short)", "Skip (too long)"])
def test_length_skipped_contig_preserves_direct_pass_and_marks_coverage_incomplete(
    contig_status,
):
    direct, assemblies, manifest = fixture_data(contig_status=contig_status)
    _, summary = combine(direct, assemblies, manifest)
    assert [row["final_status"] for row in summary["oligos"]] == [
        "Incomplete",
        "Flag",
        "Incomplete",
    ]
    assert all(row["assembly_coverage"] == "incomplete" for row in summary["oligos"])
    assert summary["flagged_oligos"] == 1


def test_deliberate_taxonomy_skip_is_complete_assembly_coverage():
    direct, assemblies, manifest = fixture_data(
        hit_status="Flag (Cleared)", contig_status="Pass (Skipped Taxonomy)"
    )
    _, summary = combine(direct, assemblies, manifest)
    assert all(row["assembly_coverage"] == "screened" for row in summary["oligos"])
    assert not any(row["incomplete_screening"] for row in summary["oligos"])


def test_excluded_and_isolated_sources_remain_direct_only_with_reasons():
    direct, assemblies, manifest = fixture_data(
        hit_status="Flag (Cleared)", contig_status="Pass"
    )
    manifest["exclusions"] = [{"id": "q0", "reason": "Non-ACGT symbol"}]
    manifest["contigs"] = []
    assemblies["queries"] = {}
    _, summary = combine(direct, assemblies, manifest)
    assert summary["oligos"][0]["assembly_coverage"] == "excluded"
    assert summary["oligos"][0]["final_status"] == "Pass"
    assert summary["oligos"][1]["assembly_coverage"] == "no_contig"
    assert summary["oligos"][1]["final_status"] == "Pass"


def test_assembly_failure_does_not_replace_direct_result():
    direct, assemblies, manifest = fixture_data(
        hit_status="Flag (Cleared)", contig_status="Pass"
    )
    manifest["contigs"] = []
    manifest["error"] = "No usable sequences"
    _, summary = combine(direct, None, manifest)
    assert [row["final_status"] for row in summary["oligos"]] == ["Pass"] * 3
    assert all(row["assembly_coverage"] == "failed" for row in summary["oligos"])
    assert all(row["incomplete_screening"] for row in summary["oligos"])


def test_nested_arguments_force_single_assembly_pass():
    from types import SimpleNamespace
    from unittest.mock import patch

    from commec.screen import Screen

    outer = Screen.__new__(Screen)
    outer.assembly_contigs = [{"id": "a0"}]
    outer.assembly_manifest = {"error": ""}
    outer.params = SimpleNamespace(directory_prefix="/tmp/example")
    outer.screen_data = SimpleNamespace(assembly_info=None)
    outer.assembly_screen = None
    args = SimpleNamespace(
        fasta_file="input",
        output_prefix="out",
        force=False,
        resume=True,
        assemble_oligos=True,
        user_specified_args={"resume"},
    )
    captured = {}
    with patch.object(
        Screen, "run", lambda self, nested: captured.update(vars(nested))
    ):
        outer.screen_oligo_assemblies(args)
    assert captured["force"] is True
    assert captured["resume"] is False
    assert captured["assemble_oligos"] is False
    assert {"force", "resume", "assemble_oligos"}.issubset(
        captured["user_specified_args"]
    )
    assert outer.assembly_manifest["contig_screen_time"] == ""


def test_nested_system_exit_restores_parent_handlers_and_blocks_export():
    import logging
    from types import SimpleNamespace
    from unittest.mock import patch

    from commec.screen import Screen

    outer = Screen.__new__(Screen)
    outer.assembly_contigs = [{"id": "a0"}]
    outer.assembly_manifest = {"error": ""}
    outer.params = SimpleNamespace(directory_prefix="/tmp/example")
    outer.screen_data = SimpleNamespace(assembly_info=None)
    outer.assembly_screen = None
    args = SimpleNamespace(
        fasta_file="input",
        output_prefix="out",
        force=False,
        resume=True,
        assemble_oligos=True,
        user_specified_args={"resume"},
    )
    logger = logging.getLogger("commec")
    parent = logging.StreamHandler()
    outer._parent_log_handlers = [parent]
    nested = logging.FileHandler("/tmp/oligo-nested-test.log")
    logger.addHandler(nested)
    with patch.object(Screen, "run", side_effect=SystemExit(1)):
        try:
            outer.screen_oligo_assemblies(args)
        except RuntimeError:
            pass
    assert parent in logger.handlers
    assert nested not in logger.handlers
    assert nested.stream is None
    logger.removeHandler(parent)


def test_parent_early_exit_without_manifest_preserves_direct_results():
    from types import SimpleNamespace

    from commec.config.result import QueryResult, ScreenResult
    from commec.screen import Screen

    screen = Screen()
    screen.params = SimpleNamespace(config={"assemble_oligos": True})
    screen.queries = {}
    screen.screen_data = ScreenResult(queries={"q": QueryResult("q")})
    screen.finalize_assembly_results()
    assert screen.screen_data.assembly_info.status == "error"
    assert screen.screen_data.queries["q"].assembly.coverage_status == "failed"


def test_early_exit_finalizes_direct_status_before_failed_assembly():
    from types import SimpleNamespace

    from commec.config.result import QueryResult, ScreenResult
    from commec.screen import Screen

    screen = Screen()
    screen.params = SimpleNamespace(config={"assemble_oligos": True})
    screen.queries = {}
    screen.screen_data = ScreenResult(queries={"q": QueryResult("q")})
    screen.finalize_assembly_results()
    assert screen.screen_data.assembly_info.status == "error"
    assert screen.screen_data.queries["q"].assembly.coverage_status == "failed"


def test_missing_contig_results_are_failed_not_no_contig():
    direct, _, manifest = fixture_data()
    _, summary = combine(direct, None, manifest)
    assert all(row["assembly_coverage"] == "failed" for row in summary["oligos"])
    assert all(row["incomplete_screening"] for row in summary["oligos"])


def test_recursive_flag_search_excludes_internal_assembly_results(tmp_path):
    import json

    from commec.flag import run

    result = {
        "commec_info": {"json_output_version": "0.6"},
        "queries": {
            "q": {"query": "q", "status": {"screen_status": "Pass"}, "hits": []}
        },
    }
    (tmp_path / "plain.output.json").write_text(json.dumps(result))
    internal = tmp_path / ".assembly_work"
    internal.mkdir()
    (internal / "contig.output.json").write_text(json.dumps(result))
    output = tmp_path / "out"
    args = type(
        "Args",
        (),
        {
            "directory": str(tmp_path),
            "recursive": True,
            "output": str(output),
            "evalportal_format": False,
        },
    )()
    run(args)
    text = (output / "screen_pipeline_status.csv").read_text()
    assert text.count("q,") == 1


def test_direct_pass_with_new_assembly_flag_and_typed_result():
    from commec.config.json_io import encode_dict_to_screen_data

    direct, assemblies, manifest = fixture_data()
    effective, summary = combine(
        encode_dict_to_screen_data(direct),
        encode_dict_to_screen_data(assemblies),
        manifest,
    )
    query = effective.queries["q1"]
    assert query.status.screen_status == "Pass"
    assert query.combined_status == "Flag"
    assert query.assembly.evidence[0].hit.region.query_start == 31
    assert query.hits == []


def test_nested_taxonomy_coordinates_are_projected_and_reverse_strand_is_correct():
    direct, assemblies, manifest = fixture_data()
    hit = assemblies["queries"]["a0"]["hits"][0]
    hit["annotations"] = {
        "controlled_taxonomy": {
            "controlled_taxa": [{"query_start": 111, "query_end": 150}],
            "non-controlled_taxa": [],
        }
    }
    manifest["contigs"][0]["placements"][1]["strand"] = "-"
    effective, _ = combine(direct, assemblies, manifest)
    evidence = effective["queries"]["q1"]["assembly"]["evidence"][0]
    taxon = evidence["hit"]["annotations"]["controlled_taxonomy"]["controlled_taxa"][0]
    assert evidence["hit"]["region"]["query_start"] == 31
    assert evidence["hit"]["region"]["query_end"] == 70
    assert taxon["query_start"] == 31
    assert taxon["query_end"] == 70


def test_empty_input_summary_has_zero_percentage():
    _, summary = combine(
        {"queries": {}}, None, {"records": [], "contigs": [], "exclusions": []}
    )
    assert summary["original_oligos"] == 0
    assert summary["flagged_oligo_percentage"] == 0.0


def test_nested_taxonomy_is_clipped_to_the_hit_intersection():
    direct, assemblies, manifest = fixture_data()
    hit = assemblies["queries"]["a0"]["hits"][0]
    hit["annotations"] = {
        "controlled_taxonomy": {
            "controlled_taxa": [{"query_start": 111, "query_end": 150}],
            "non-controlled_taxa": [],
        }
    }
    manifest["contigs"][0]["placements"][1] = {
        "id": "q1",
        "strand": "-",
        "start": 80,
        "end": 180,
    }
    effective, _ = combine(direct, assemblies, manifest)
    evidence = effective["queries"]["q1"]["assembly"]["evidence"][0]
    assert evidence["hit"]["region"]["query_start"] == 31
    assert evidence["hit"]["region"]["query_end"] == 70
    taxon = evidence["hit"]["annotations"]["controlled_taxonomy"]["controlled_taxa"][0]
    assert taxon["query_start"] == 31
    assert taxon["query_end"] == 70


def test_console_summary_uses_combined_status():
    from commec.config.json_io import encode_dict_to_screen_data

    direct, assemblies, manifest = fixture_data()
    result, _ = combine(
        encode_dict_to_screen_data(direct),
        encode_dict_to_screen_data(assemblies),
        manifest,
    )
    assert result.get_flag_data()["overall"].tolist() == ["Pass", "Flag", "Pass"]
    assert result.get_rationale_data()["overall"].tolist() == ["Pass", "Flag", "Pass"]


def test_combined_status_consumers_use_assembly_evidence(tmp_path):
    from commec.flag import read_flags_from_json
    from commec.gui.server import _summarize_output
    from commec.oligo_assembly import write_json

    direct, assemblies, manifest = fixture_data()
    direct["queries"]["q1"]["status"].update(
        {
            "screen_status": "Pass",
            "biorisk": "Pass",
            "protein_taxonomy": "Pass",
            "nucleotide_taxonomy": "Pass",
            "low_concern": "Skip",
        }
    )
    hit = assemblies["queries"]["a0"]["hits"][0]
    hit["annotations"] = {
        "controlled_taxonomy": {
            "statistics": {
                "controlled_bacteria": 0,
                "controlled_viruses": 1,
                "controlled_fungi": 0,
                "controlled_parasites": 0,
                "number_of_non-controlled_taxids": 0,
            }
        }
    }
    effective, _ = combine(direct, assemblies, manifest)
    effective["commec_info"] = {"json_output_version": "0.7"}
    write_json(tmp_path / "combined.output.json", effective)
    flagged = read_flags_from_json(tmp_path / "combined.output.json")
    assert flagged[1]["flag"] == "Flag"
    assert flagged[1]["virus_flag"] is True
    summary = _summarize_output(tmp_path)
    assert summary["queries"][1]["status"] == "Flag"
    assert "Assembly coverage: screened" in summary["queries"][1]["rationale"]


def test_disabled_json_keeps_legacy_shape(tmp_path):
    from commec.config.json_io import (
        encode_screen_data_to_json,
        get_screen_data_from_json,
    )
    from commec.config.result import ScreenResult

    path = tmp_path / "plain.output.json"
    encode_screen_data_to_json(ScreenResult(), path)
    document = path.read_text()
    assert '"assembly_info"' not in document
    assert '"assembly_summary"' not in document
    assert '"assembly"' not in document
    assert '"combined_status"' not in document
    loaded = get_screen_data_from_json(path)
    assert loaded.assembly_info is None
    assert loaded.assembly_summary is None
