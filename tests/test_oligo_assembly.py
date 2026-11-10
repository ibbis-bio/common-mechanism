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
        "Incomplete",
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
    assert [row["final_status"] for row in summary["oligos"]] == ["Incomplete"] * 3
    assert all(row["assembly_coverage"] == "incomplete" for row in summary["oligos"])
    assert summary["flagged_oligos"] == 0


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
