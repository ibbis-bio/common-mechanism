"""Assemble oligos and attribute contig concerns to their source records."""

import argparse
import copy
import csv
import hashlib
import json
import math
import os
import subprocess
import sys
from collections import Counter
from dataclasses import asdict, is_dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from commec.config.json_io import encode_dict_to_screen_data
from commec.config.result import (
    ASSEMBLY_JSON_FORMAT_VERSION,
    AssemblyEvidence,
    AssemblyQueryResult,
    AssemblyRunInfo,
    ScreenStatus,
    ScreenStep,
)

ACTIVE = {ScreenStatus.FLAG, ScreenStatus.WARN}
INCOMPLETE = {ScreenStatus.ERROR, ScreenStatus.STOP, ScreenStatus.NULL}
COVERAGE_INCOMPLETE = INCOMPLETE | {
    ScreenStatus.SKIP,
    ScreenStatus.SKIP_SHORT,
    ScreenStatus.SKIP_LONG,
}


def write_json(path, value):
    def clean(item):
        if isinstance(item, Counter):
            return dict(item)
        if is_dataclass(item):
            return clean(asdict(item))
        if isinstance(item, float) and not math.isfinite(item):
            return None
        if isinstance(item, dict):
            return {key: clean(val) for key, val in item.items()}
        if isinstance(item, list):
            return [clean(val) for val in item]
        return item

    Path(path).write_text(
        json.dumps(clean(value), indent=2, allow_nan=False, default=str) + "\n"
    )


def load_records(path):
    """Keep record identity separate from possibly repeated FASTA headers."""
    records = []
    for index, record in enumerate(SeqIO.parse(path, "fasta")):
        records.append(
            {
                "id": f"q{index:06d}",
                "header": record.description,
                "original_id": record.id,
                "sequence": str(record.seq).upper(),
                "record_index": index,
            }
        )
    if not records:
        raise ValueError("The input contains no FASTA records.")
    return records


def write_fasta(path, records):
    Path(path).write_text(
        "".join(f">{record['id']}\n{record['sequence']}\n" for record in records)
    )


def assemble(records, min_overlap=15):
    """Replay each actual walk and verify every source placement."""
    import oligraph

    valid, exclusions = {}, []
    for record in records:
        sequence = record["sequence"]
        reason = None
        if not sequence or set(sequence) - set("ACGT"):
            reason = "Empty sequence or non-ACGT symbol"
        elif len(sequence) > oligraph.MAX_SEQUENCE_LENGTH:
            reason = "Sequence exceeds the assembler input length limit"
        if reason:
            exclusions.append({"id": record["id"], "reason": reason})
        else:
            valid[record["id"]] = sequence
    if not valid:
        return {
            "assembler": "oligraph",
            "version": oligraph.__version__,
            "min_overlap": min_overlap,
            "matching": "exact",
            "method": "all",
            "coordinate_system": "zero-based half-open source placements",
            "contigs": [],
            "exclusions": exclusions,
            "stats": {
                "n_sequences": 0,
                "n_skipped": len(exclusions),
                "len_min": 0,
                "len_max": 0,
                "n_edges": 0,
                "n_isolated": 0,
                "n_contigs": 0,
            },
            "records": [
                {k: v for k, v in record.items() if k != "sequence"}
                | {"length": len(record["sequence"])}
                for record in records
            ],
        }, []
    graph = oligraph.build(
        valid, min_overlap=min_overlap, method="all", on_invalid="raise"
    )
    edge_map = {}
    for edge in graph.edges:
        a, b = edge.from_name, edge.to_name
        sa, sb = edge.from_strand.symbol, edge.to_strand.symbol
        edge_map[a, sa, b, sb] = edge.overlap_len
        edge_map[b, flip(sb), a, flip(sa)] = edge.overlap_len
    contigs, contig_records = [], []
    for index, contig in enumerate(graph.contigs):
        assembled, placements = "", []
        previous = None
        for identifier, strand in contig.path:
            symbol = strand.symbol
            sequence = valid[identifier]
            oriented = (
                sequence if symbol == "+" else str(Seq(sequence).reverse_complement())
            )
            overlap = (
                0 if previous is None else edge_map[(*previous, identifier, symbol)]
            )
            if overlap and assembled[-overlap:] != oriented[:overlap]:
                raise ValueError(
                    "The reported assembly edge does not match its sequences."
                )
            start = len(assembled) - overlap
            assembled += oriented[overlap:]
            placements.append(
                {
                    "id": identifier,
                    "strand": symbol,
                    "start": start,
                    "end": start + len(oriented),
                }
            )
            previous = identifier, symbol
        if assembled != contig.sequence:
            raise ValueError("The source placements do not reproduce the contig.")
        for placement in placements:
            sequence = valid[placement["id"]]
            if placement["strand"] == "-":
                sequence = str(Seq(sequence).reverse_complement())
            if assembled[placement["start"] : placement["end"]] != sequence:
                raise ValueError("A source oligo does not match its contig placement.")
        identifier = f"a{index:06d}"
        contig_records.append({"id": identifier, "sequence": assembled})
        contigs.append(
            {
                "id": identifier,
                "length": len(assembled),
                "sha256": hashlib.sha256(assembled.encode()).hexdigest(),
                "component": contig.component,
                "topology": str(contig.topology),
                "branches": contig.branches,
                "placements": placements,
            }
        )
    manifest = {
        "assembler": "oligraph",
        "version": oligraph.__version__,
        "min_overlap": min_overlap,
        "matching": "exact",
        "method": "all",
        "coordinate_system": "zero-based half-open source placements",
        "contigs": contigs,
        "exclusions": exclusions,
        "stats": graph.stats.to_dict(),
        "records": [
            {k: v for k, v in r.items() if k != "sequence"}
            | {"length": len(r["sequence"])}
            for r in records
        ],
    }
    return manifest, contig_records


def flip(strand):
    return "-" if strand == "+" else "+"


def project_region(placement, start, end):
    """Return the original oligo interval intersecting a contig region."""
    left, right = max(start, placement["start"]), min(end, placement["end"])
    if left >= right:
        return None
    if placement["strand"] == "+":
        return left - placement["start"], right - placement["start"]
    return placement["end"] - right, placement["end"] - left


def _assembly_evidence(contig, placement, hit, start, end):
    """Project a contig hit and its nested match coordinates to an oligo."""
    original_start, original_end = project_region(placement, start, end)
    projected = copy.deepcopy(hit)
    projected.region.query_start = original_start + 1
    projected.region.query_end = original_end
    taxonomy = projected.annotations.get("controlled_taxonomy", {})
    for group in ("controlled_taxa", "non-controlled_taxa"):
        retained = []
        for taxon in taxonomy.get(group, []):
            left = taxon.get("query_start")
            right = taxon.get("query_end")
            if not all(
                isinstance(x, int) and not isinstance(x, bool) for x in (left, right)
            ):
                continue
            region = project_region(placement, min(left, right) - 1, max(left, right))
            if region is None:
                continue
            taxon["query_start"], taxon["query_end"] = region[0] + 1, region[1]
            retained.append(taxon)
        taxonomy[group] = retained
    projected.annotations["assembly_origin"] = {
        "contig": contig["id"],
        "contig_region": [start + 1, end],
        "evidence_scope": "Contig-derived evidence, not an independent oligo alignment",
    }
    return AssemblyEvidence(
        contig_id=contig["id"],
        contig_region=[start + 1, end],
        oligo_region=[original_start + 1, original_end],
        strand=placement["strand"],
        status=hit.recommendation.status,
        hit=projected,
    )


def _initial_coverage(manifest, assemblies):
    """Record which original records occur in each returned contig."""
    excluded = {item["id"]: item["reason"] for item in manifest["exclusions"]}
    failure = manifest.get("error")
    coverage = {}
    for record in manifest["records"]:
        identifier = record["id"]
        status = "excluded" if identifier in excluded else "no_contig"
        reason = excluded.get(
            identifier, "No returned contig contains this original oligo."
        )
        if failure and identifier not in excluded:
            status, reason = "failed", failure
        coverage[identifier] = AssemblyQueryResult(
            enabled=True, coverage_status=status, coverage_reason=reason
        )
    if assemblies is None:
        return coverage
    for contig in manifest["contigs"]:
        result = assemblies.queries[contig["id"]]
        status = result.status.screen_status
        for placement in contig["placements"]:
            query = coverage[placement["id"]]
            query.contig_ids.append(contig["id"])
            if status in COVERAGE_INCOMPLETE:
                query.coverage_status = "incomplete"
                query.coverage_reason = f"Contig {contig['id']} screening: {status}."
            elif query.coverage_status != "incomplete":
                query.coverage_status = "screened"
                query.coverage_reason = ""
    return coverage


def _combined_status(direct_status, assembly):
    """Keep known biological evidence separate from incomplete screening."""
    inherited = {item.status for item in assembly.evidence}
    if direct_status == ScreenStatus.FLAG or ScreenStatus.FLAG in inherited:
        return ScreenStatus.FLAG
    if direct_status == ScreenStatus.WARN or ScreenStatus.WARN in inherited:
        return ScreenStatus.WARN
    if direct_status == ScreenStatus.ERROR:
        return direct_status
    if assembly.coverage_status == "incomplete":
        return ScreenStatus.STOP
    return direct_status


def _restore_contig_coverage(coverage, assemblies, contigs):
    """Restore screened coverage after active-hit attribution."""
    if assemblies is None:
        return
    for contig in contigs:
        if assemblies.queries[contig["id"]].status.screen_status in COVERAGE_INCOMPLETE:
            continue
        for placement in contig["placements"]:
            query = coverage[placement["id"]]
            if query.coverage_status == "screened":
                query.coverage_reason = (
                    "Only uncleared contig-hit intersections create assembly evidence."
                )


def combine(direct, assemblies, manifest):
    """Combine two ScreenResults, or two legacy dictionaries, without losing evidence."""
    legacy = isinstance(direct, dict)
    if legacy:
        for document in (direct, assemblies):
            if document is None:
                continue
            for query in document.get("queries", {}).values():
                for hit in query.get("hits", []):
                    hit.setdefault("name", "")
                    hit.setdefault("description", "")
                    hit.setdefault("region", {})
                    hit.setdefault("annotations", {})
                    hit.setdefault("recommendation", {}).setdefault("status", "-")
                    hit["recommendation"].setdefault("from_step", ScreenStep.BIORISK)
        direct = encode_dict_to_screen_data(direct)
        assemblies = (
            encode_dict_to_screen_data(assemblies) if assemblies is not None else None
        )
    effective = copy.deepcopy(direct)
    effective.commec_info.json_output_version = ASSEMBLY_JSON_FORMAT_VERSION
    identifiers = {record["id"] for record in manifest["records"]}
    if set(direct.queries) != identifiers:
        raise ValueError("Direct screening lost or changed an input record identity.")
    if assemblies is not None and set(assemblies.queries) != {
        c["id"] for c in manifest["contigs"]
    }:
        raise ValueError("Assembly screening did not account for every contig.")
    coverage = _initial_coverage(manifest, assemblies)
    unresolved, contig_statuses = [], {}
    flagged_contigs = set()
    for contig in manifest["contigs"]:
        if assemblies is None:
            break
        result = assemblies.queries[contig["id"]]
        final = result.status.screen_status
        contig_statuses[contig["id"]] = str(final)
        if final == ScreenStatus.FLAG:
            flagged_contigs.add(contig["id"])
        if final not in ACTIVE | COVERAGE_INCOMPLETE:
            continue
        active_hits = [h for h in result.hits if h.recommendation.status in ACTIVE]
        if not active_hits and final in ACTIVE:
            unresolved.append(
                {"contig": contig["id"], "reason": "No active hit coordinates"}
            )
        for hit in active_hits:
            if hit.recommendation.status == ScreenStatus.FLAG:
                flagged_contigs.add(contig["id"])
            left, right = hit.region.query_start, hit.region.query_end
            if not all(
                isinstance(x, int) and not isinstance(x, bool) for x in (left, right)
            ):
                unresolved.append(
                    {"contig": contig["id"], "reason": "Missing hit coordinates"}
                )
                continue
            start, end = min(left, right) - 1, max(left, right)
            if start < 0 or end > contig["length"] or start >= end:
                unresolved.append(
                    {"contig": contig["id"], "reason": "Hit coordinates out of range"}
                )
                continue
            for placement in contig["placements"]:
                if project_region(placement, start, end) is not None:
                    coverage[placement["id"]].evidence.append(
                        _assembly_evidence(contig, placement, hit, start, end)
                    )
    for item in unresolved:
        contig = next(c for c in manifest["contigs"] if c["id"] == item["contig"])
        for placement in contig["placements"]:
            coverage[placement["id"]].coverage_status = "incomplete"
            coverage[placement["id"]].coverage_reason = item["reason"]
    _restore_contig_coverage(coverage, assemblies, manifest["contigs"])
    rows = []
    for record in manifest["records"]:
        identifier = record["id"]
        assembly = coverage[identifier]
        query = effective.queries[identifier]
        query.assembly = assembly
        query.combined_status = _combined_status(query.status.screen_status, assembly)
        query.description = record["header"]
        inherited = {item.status for item in assembly.evidence}
        incomplete = (
            query.status.screen_status in INCOMPLETE
            or assembly.coverage_status in {"incomplete", "failed"}
        )
        rows.append(
            {
                "query_id": identifier,
                "original_id": record["original_id"],
                "record_index": record["record_index"],
                "direct_status": str(query.status.screen_status),
                "assembly_status": "Flag"
                if ScreenStatus.FLAG in inherited
                else "Warning"
                if inherited
                else "No attributed hit",
                "assembly_coverage": assembly.coverage_status,
                "final_status": str(query.combined_status),
                "incomplete_screening": incomplete,
                "contigs": sorted(assembly.contig_ids),
                "evidence": [item.to_summary_dict() for item in assembly.evidence],
            }
        )
    effective.assembly_info = AssemblyRunInfo(
        enabled=True,
        status="error"
        if manifest.get("error")
        else "incomplete"
        if unresolved or any(r["assembly_coverage"] == "incomplete" for r in rows)
        else "complete",
        assembler=manifest.get("assembler", "oligraph"),
        version=manifest.get("version", ""),
        min_overlap=manifest.get("min_overlap", 15),
        matching=manifest.get("matching", "exact"),
        method=manifest.get("method", "all"),
        contigs=[
            {**c, "status": contig_statuses.get(c["id"], "Not screened")}
            for c in manifest["contigs"]
        ],
        exclusions=manifest["exclusions"],
        unresolved_attributions=unresolved,
        error=manifest.get("error", ""),
        contig_screen_time=manifest.get("contig_screen_time", ""),
        manifest=manifest.get("path", ""),
    )
    for query in effective.queries.values():
        if not query.assembly.contig_ids:
            continue
        query.assembly.contig_ids = sorted(set(query.assembly.contig_ids))
    effective.assembly_summary = {
        "original_oligos": len(rows),
        "assemblies": len(manifest["contigs"]),
        "flagged_assemblies": len(flagged_contigs),
        "direct_status_counts": dict(Counter(row["direct_status"] for row in rows)),
        "final_status_counts": dict(Counter(row["final_status"] for row in rows)),
        "assembly_coverage_counts": dict(
            Counter(row["assembly_coverage"] for row in rows)
        ),
        "flagged_oligos": sum(
            row["final_status"] == str(ScreenStatus.FLAG) for row in rows
        ),
        "newly_flagged_oligos": sum(
            ScreenStatus.FLAG
            in {item.status for item in coverage[row["query_id"]].evidence}
            and not any(
                h.recommendation.status == ScreenStatus.FLAG
                for h in direct.queries[row["query_id"]].hits
            )
            for row in rows
        ),
        "unresolved_attributions": unresolved,
        "assembly_status_counts": dict(Counter(contig_statuses.values())),
        "assembly_statuses": contig_statuses,
        "assembly_exclusions": manifest["exclusions"],
        "oligos": rows,
        "pool_verdict": None,
        "attribution_rule": "Nonempty intersection with an uncleared final hit in an actual contig path",
    }
    effective.assembly_summary["flagged_oligo_percentage"] = (
        100 * effective.assembly_summary["flagged_oligos"] / len(rows)
    )
    return (asdict(effective) if legacy else effective), effective.assembly_summary


def screen(fasta, directory, databases, threads):
    directory.mkdir(parents=True, exist_ok=True)
    command = [
        sys.executable,
        "-m",
        "commec.cli",
        "screen",
        "-d",
        str(databases),
        "-t",
        str(threads),
        "--force",
        "-o",
        str(directory),
        str(fasta),
    ]
    source_root = str(Path(__file__).resolve().parent.parent)
    env = dict(os.environ)
    env["PYTHONPATH"] = source_root + (
        os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else ""
    )
    result = subprocess.run(
        command, capture_output=True, text=True, timeout=7200, env=env
    )
    (directory / "process.stdout.txt").write_text(result.stdout)
    (directory / "process.stderr.txt").write_text(result.stderr)
    if result.returncode:
        raise RuntimeError(
            f"COMMEC failed with exit code {result.returncode}. See process logs."
        )
    candidates = []
    for path in directory.glob("*.json"):
        try:
            document = json.loads(path.read_text())
        except (ValueError, OSError):
            continue
        if "queries" in document and "commec_info" in document:
            candidates.append((path, document))
    if len(candidates) != 1:
        raise RuntimeError("Expected exactly one COMMEC result JSON.")
    return candidates[0][1]


def render(effective_path, summary, manifest, output):
    """Render assembly evidence through the normal interactive report."""
    from commec.config.json_io import get_screen_data_from_json
    from commec.utils.json_html_output import build_report_model, render_report_model

    model = build_report_model(get_screen_data_from_json(str(effective_path)))
    names = {
        record["id"]: f"{record['original_id']} (record {record['record_index'] + 1})"
        for record in manifest["records"]
        if "original_id" in record
    }
    for sequence in model["sequences"]:
        sequence["name"] = names.get(sequence["id"], sequence["name"])
    Path(output).write_text(render_report_model(model))


def run(input_path, output_dir, databases, threads=8):
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    records = load_records(input_path)
    original_fasta = output / "originals.fasta"
    write_fasta(original_fasta, records)
    direct = screen(original_fasta, output / "direct", databases, threads)
    try:
        manifest, contigs = assemble(records)
    except Exception as error:
        manifest = {
            "assembler": "oligraph",
            "version": "",
            "min_overlap": 15,
            "matching": "exact",
            "method": "all",
            "records": records,
            "contigs": [],
            "exclusions": [],
            "error": str(error),
        }
        contigs = []
    write_json(output / "assembly_manifest.json", manifest)
    contig_fasta = output / "assemblies.fasta"
    write_fasta(contig_fasta, contigs)
    assembled = (
        screen(contig_fasta, output / "contigs", databases, threads)
        if contigs
        else None
    )
    effective, summary = combine(direct, assembled if contigs else None, manifest)
    write_json(output / "effective_report.json", effective)
    write_json(output / "summary.json", summary)
    with (output / "oligos.csv").open("w", newline="") as stream:
        fields = [
            "query_id",
            "original_id",
            "record_index",
            "direct_status",
            "assembly_status",
            "assembly_coverage",
            "final_status",
            "incomplete_screening",
            "contigs",
        ]
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(summary["oligos"])
    render(output / "effective_report.json", summary, manifest, output / "report.html")
    return {
        key: value
        for key, value in summary.items()
        if key not in {"oligos", "assembly_statuses"}
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--databases", required=True)
    parser.add_argument("--threads", type=int, default=8)
    arguments = parser.parse_args()
    print(
        json.dumps(
            run(
                arguments.input,
                arguments.output,
                arguments.databases,
                arguments.threads,
            )
        )
    )
