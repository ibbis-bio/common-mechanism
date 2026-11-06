"""Verify source attribution without using private screening sequences."""

import copy

import pytest

from commec.oligo_assembly_poc import assemble, combine, load_records, project_region


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
    assert effective["queries"]["q1"]["hits"][0]["region"] == {
        "query_start": 31,
        "query_end": 70,
    }


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
    assert all(r["incomplete_screening"] for r in summary["oligos"])
    assert summary["flagged_assemblies"] == 1
    assert summary["assembly_statuses"] == {"a0": "Error"}


def test_direct_error_keeps_known_flag_evidence():
    direct, assemblies, manifest = fixture_data(
        hit_status="Flag (Cleared)", contig_status="Pass"
    )
    direct["queries"]["q0"]["status"]["screen_status"] = "Error"
    direct["queries"]["q0"]["hits"] = [{"recommendation": {"status": "Flag"}}]
    _, summary = combine(direct, assemblies, manifest)
    assert summary["oligos"][0]["final_status"] == "Flag"
    assert summary["oligos"][0]["incomplete_screening"]


def test_html_shows_duplicate_record_numbers_and_assembly_skip(tmp_path):
    from dataclasses import asdict

    from commec.config.result import ScreenResult
    from commec.oligo_assembly_poc import render, write_json

    direct, assemblies, manifest = fixture_data(
        hit_status="Flag (Cleared)", contig_status="Skip (too short)"
    )
    effective, summary = combine(direct, assemblies, manifest)
    effective_path, output = tmp_path / "effective.json", tmp_path / "report.html"
    document = asdict(ScreenResult())
    document["queries"] = effective["queries"]
    write_json(effective_path, document)
    render(effective_path, summary, manifest, output)
    page = output.read_text()
    assert "same_header (record 1)" in page
    assert "same_header (record 2)" in page
    assert "Every assembly and its screening status" in page
    assert "Skip (too short)" in page
