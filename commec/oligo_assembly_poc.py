"""Trace assembly screening evidence to original oligos in a small proof of concept."""

import argparse
import copy
import csv
import hashlib
import html
import json
import math
import subprocess
import sys
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

ACTIVE = {"Flag", "Warning"}
INCOMPLETE = {"Error", "Incomplete", "-"}


def write_json(path, value):
    def clean(item):
        if isinstance(item, float) and not math.isfinite(item):
            return None
        if isinstance(item, dict):
            return {key: clean(val) for key, val in item.items()}
        if isinstance(item, list):
            return [clean(val) for val in item]
        return item

    Path(path).write_text(json.dumps(clean(value), indent=2, allow_nan=False) + "\n")


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


def combine(direct, assemblies, manifest):
    """Use final uncleared hits, keeping raw direct results unchanged."""
    effective = copy.deepcopy(direct)
    identifiers = {record["id"] for record in manifest["records"]}
    if set(direct["queries"]) != identifiers:
        raise ValueError("Direct screening lost or changed an input record identity.")
    if set(assemblies["queries"]) != {c["id"] for c in manifest["contigs"]}:
        raise ValueError("Assembly screening did not account for every contig.")
    evidence = {identifier: [] for identifier in identifiers}
    unresolved, contig_statuses = [], {}
    incomplete_sources = set()
    flagged_contigs = set()
    for contig in manifest["contigs"]:
        result = assemblies["queries"][contig["id"]]
        final = result["status"]["screen_status"]
        contig_statuses[contig["id"]] = final
        if final in INCOMPLETE:
            incomplete_sources.update(p["id"] for p in contig["placements"])
        if final == "Flag":
            flagged_contigs.add(contig["id"])
        if final not in ACTIVE | INCOMPLETE:
            continue
        active_hits = [
            h for h in result.get("hits", []) if h["recommendation"]["status"] in ACTIVE
        ]
        if not active_hits and final in ACTIVE:
            incomplete_sources.update(p["id"] for p in contig["placements"])
            unresolved.append(
                {"contig": contig["id"], "reason": "No active hit coordinates"}
            )
        for hit in active_hits:
            if hit["recommendation"]["status"] == "Flag":
                flagged_contigs.add(contig["id"])
            region = hit.get("region", {})
            bounds = region.get("query_start"), region.get("query_end")
            if not all(isinstance(x, int) and not isinstance(x, bool) for x in bounds):
                incomplete_sources.update(p["id"] for p in contig["placements"])
                unresolved.append(
                    {"contig": contig["id"], "reason": "Missing hit coordinates"}
                )
                continue
            start, end = min(bounds) - 1, max(bounds)
            if start < 0 or end > contig["length"] or start >= end:
                incomplete_sources.update(p["id"] for p in contig["placements"])
                unresolved.append(
                    {"contig": contig["id"], "reason": "Hit coordinates out of range"}
                )
                continue
            for placement in contig["placements"]:
                projected = project_region(placement, start, end)
                if projected is None:
                    continue
                original_start, original_end = projected
                trace = {
                    "contig": contig["id"],
                    "contig_region": [start, end],
                    "oligo_region": [original_start, original_end],
                    "strand": placement["strand"],
                    "status": hit["recommendation"]["status"],
                    "hit": hit["name"],
                }
                evidence[placement["id"]].append(trace)
                projected_hit = copy.deepcopy(hit)
                projected_hit["region"]["query_start"] = original_start + 1
                projected_hit["region"]["query_end"] = original_end
                projected_hit["name"] = f"{hit['name']} [assembly {contig['id']}]"
                projected_hit.setdefault("annotations", {})["assembly_origin"] = trace
                effective["queries"][placement["id"]]["hits"].append(projected_hit)
    rows = []
    for record in manifest["records"]:
        identifier = record["id"]
        query = effective["queries"][identifier]
        direct_status = direct["queries"][identifier]["status"]["screen_status"]
        direct_active = {
            hit["recommendation"]["status"]
            for hit in direct["queries"][identifier].get("hits", [])
            if direct_status in INCOMPLETE and hit["recommendation"]["status"] in ACTIVE
        }
        inherited = {item["status"] for item in evidence[identifier]}
        final = direct_status
        if direct_status == "Flag" or "Flag" in inherited | direct_active:
            final = "Flag"
        elif direct_status == "Warning" or "Warning" in inherited | direct_active:
            final = "Warning"
        elif identifier in incomplete_sources and direct_status not in INCOMPLETE:
            final = "Incomplete"
        query["status"]["screen_status"] = final
        if evidence[identifier]:
            ids = sorted({item["contig"] for item in evidence[identifier]})
            query["status"]["rationale"] = (
                f"Assembly evidence from {', '.join(ids)}. Direct screening: {direct_status}."
            )
        if identifier in incomplete_sources:
            query["status"]["rationale"] = (
                str(query["status"].get("rationale", ""))
                + " Assembly screening was incomplete."
            ).strip()
        query["description"] = record["header"]
        rows.append(
            {
                "query_id": identifier,
                "original_id": record["original_id"],
                "record_index": record["record_index"],
                "direct_status": direct_status,
                "assembly_status": "Flag"
                if "Flag" in inherited
                else "Warning"
                if inherited
                else "No attributed hit",
                "final_status": final,
                "incomplete_screening": identifier in incomplete_sources
                or direct_status in INCOMPLETE,
                "contigs": sorted({item["contig"] for item in evidence[identifier]}),
                "evidence": evidence[identifier],
            }
        )
    summary = {
        "original_oligos": len(rows),
        "assemblies": len(manifest["contigs"]),
        "flagged_assemblies": len(flagged_contigs),
        "direct_status_counts": dict(Counter(row["direct_status"] for row in rows)),
        "final_status_counts": dict(Counter(row["final_status"] for row in rows)),
        "flagged_oligos": sum(row["final_status"] == "Flag" for row in rows),
        "newly_flagged_oligos": sum(
            row["final_status"] == "Flag" and row["direct_status"] != "Flag"
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
    summary["flagged_oligo_percentage"] = 100 * summary["flagged_oligos"] / len(rows)
    return effective, summary


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
    result = subprocess.run(command, capture_output=True, text=True, timeout=7200)
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
    from commec.config.json_io import get_screen_data_from_json
    from commec.utils.json_html_output import build_report_model, render_report_model

    model = build_report_model(get_screen_data_from_json(str(effective_path)))

    def label(record):
        return f"{record['original_id']} (record {record['record_index'] + 1})"

    names = {record["id"]: label(record) for record in manifest["records"]}
    for sequence in model["sequences"]:
        sequence["name"] = names[sequence["id"]]
    page = render_report_model(model)
    escaped = html.escape
    rows = "".join(
        "<tr>"
        + "".join(
            f"<td>{escaped(str(value))}</td>"
            for value in [
                label(r),
                r["direct_status"],
                r["assembly_status"],
                r["final_status"],
                ", ".join(r["contigs"]),
            ]
        )
        + "</tr>"
        for r in summary["oligos"]
    )
    traces = "".join(
        f"<li>{escaped(label(row))}: contig {escaped(item['contig'])}, "
        f"contig bases {item['contig_region'][0] + 1}-{item['contig_region'][1]}, "
        f"oligo bases {item['oligo_region'][0] + 1}-{item['oligo_region'][1]}, strand {item['strand']}.</li>"
        for row in summary["oligos"]
        for item in row["evidence"]
    )
    assembly_rows = "".join(
        f"<tr><td>{escaped(contig['id'])}</td><td>{contig['length']}</td>"
        f"<td>{escaped(summary['assembly_statuses'][contig['id']])}</td></tr>"
        for contig in manifest["contigs"]
    )
    panel = (
        '<section style="margin:24px auto;padding:20px;max-width:1240px;font-family:Arial,sans-serif;background:#f1f5fa;color:#18294a">'
        "<h2>Assembly screening proof of concept</h2>"
        f"<p>{summary['flagged_assemblies']} of {summary['assemblies']} assemblies flagged. "
        f"{summary['flagged_oligos']} of {summary['original_oligos']} original oligos flagged "
        f"({summary['flagged_oligo_percentage']:.1f} percent). No pool verdict is assigned.</p>"
        "<p>Assembly attribution identifies oligos overlapping the final region of concern. "
        "It does not assert that an oligo is independently a sequence of concern. "
        "Direct screening results are retained separately.</p>"
        '<table cellpadding="6"><thead><tr><th>Oligo</th><th>Direct</th><th>Assembly evidence</th><th>Final</th><th>Contigs</th></tr></thead>'
        f"<tbody>{rows}</tbody></table><details><summary>Attribution coordinates</summary><ul>{traces}</ul></details>"
        "<h3>Every assembly and its screening status</h3>"
        '<table cellpadding="6"><thead><tr><th>Assembly</th><th>Length</th><th>Screening status</th></tr></thead>'
        f"<tbody>{assembly_rows}</tbody></table>"
        f"<p>Unresolved attributions: {len(summary['unresolved_attributions'])}. "
        f"Excluded from assembly: {len(manifest['exclusions'])}. "
        "Skipped screening remains visible in the report. Coordinates shown here start at one.</p></section>"
    )
    Path(output).write_text(page.replace("<body>", "<body>\n" + panel, 1))


def run(input_path, output_dir, databases, threads=8):
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    records = load_records(input_path)
    original_fasta = output / "originals.fasta"
    write_fasta(original_fasta, records)
    manifest, contigs = assemble(records)
    write_json(output / "assembly_manifest.json", manifest)
    contig_fasta = output / "assemblies.fasta"
    write_fasta(contig_fasta, contigs)
    direct = screen(original_fasta, output / "direct", databases, threads)
    assembled = (
        screen(contig_fasta, output / "contigs", databases, threads)
        if contigs
        else {"queries": {}}
    )
    effective, summary = combine(direct, assembled, manifest)
    write_json(output / "effective_report.json", effective)
    write_json(output / "summary.json", summary)
    with (output / "oligos.csv").open("w", newline="") as stream:
        fields = [
            "query_id",
            "original_id",
            "record_index",
            "direct_status",
            "assembly_status",
            "final_status",
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
