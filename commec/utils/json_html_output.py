"""
Convert a Commec Screen JSON object (or a ScreenResult object) into a self-contained,
interactive HTML report ("Commec Screening Report").

Division of responsibilities (Python vs. browser)
-------------------------------------------------
This module is deliberately thin and only does two things:

  1. :func:`build_report_model` flattens the nested ``ScreenResult``
     (queries -> hits -> free-form ``annotations``) into a plain per-hit data model.
  2. :func:`render_report_html` assembles a single HTML file: it embeds that model as
     ``window.COMMEC_REPORT`` and inlines the CSS / JS / logo assets.

Almost all of the actual report *construction* happens in the browser, in
``report_assets/report.js`` (styled by ``report.css``): grouping hits per organism,
deciding each hit's disposition, the status tabs, the sequence table, the expandable
detail view (swimlane lanes, hit rail, best-target card, control-list drill-down) and
all interactivity. If you are changing what the report looks like or does, that is
almost certainly a report.js / report.css change -- this module only governs what
*data* reaches the page.

The report is a single HTML file with no external framework and no CDN dependency
(fonts are loaded from Google Fonts but degrade gracefully offline). Hit disposition
and tab status are derived (in report.js) from commec's own authoritative per-hit
status, so the report always agrees with commec's flag counts.
"""

import os
import re
import math
import json
import base64
import argparse
import logging
import importlib.resources
from urllib.parse import unquote

from commec.config.json_io import get_screen_data_from_json
from commec.config.result import ScreenResult, ScreenStep

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Small value helpers
# ---------------------------------------------------------------------------
def _is_nan(v) -> bool:
    return isinstance(v, float) and math.isnan(v)


def _clean(v):
    """Normalise pandas/JSON ``NaN`` and empty values to ``None``."""
    if v is None or _is_nan(v):
        return None
    return v


def _first(v):
    """Several annotation fields (regulated / domain / category) are one-item lists."""
    if isinstance(v, (list, tuple)):
        return v[0] if v else None
    return v


_HTML_TAG_RE = re.compile(r"<[^>]*>")


def _clean_text(v):
    """
    Tidy a display string coming from upstream databases. Some names carry markup
    (sometimes URL-encoded, e.g. ``%3C/i%3E`` for ``</i>``), stray wrapping quotes,
    or runs of whitespace. Decode, strip tags, and trim so it presents cleanly.
    """
    v = _clean(v)
    if not isinstance(v, str):
        return v
    v = unquote(v)  # %3C/i%3E -> </i>
    v = _HTML_TAG_RE.sub("", v)  # drop <i>...</i> and similar markup
    v = v.replace(" ", " ")  # non-breaking space -> space
    v = re.sub(r"\s+", " ", v).strip()  # collapse internal whitespace
    v = v.strip('"“”').strip()  # strip wrapping straight/smart quotes
    return v or None


def _fmt_evalue(e) -> str:
    """Format an e-value the way the report shows it (e.g. ``4.66e-19``, ``0``)."""
    e = _clean(e)
    if e is None:
        return ""
    try:
        e = float(e)
    except (TypeError, ValueError):
        return str(e)
    if math.isnan(e):
        return ""
    if e == 0:
        return "0"
    return f"{e:g}"


def _fmt_pct(p) -> str:
    p = _clean(p)
    if p is None:
        return None
    try:
        return f"{float(p):.1f}%"
    except (TypeError, ValueError):
        return str(p)


def _tool_short(tool_info: str) -> str:
    """
    Reduce a verbose tool banner to a short ``Name Version`` label.

    Handles the shapes commec records, e.g.::

        "# HMMER 3.4 (Aug 2023); http://hmmer.org/"      -> "HMMER 3.4"
        "blastx: 2.17.0+\\n Package: blast 2.17.0, ..."   -> "BLASTx 2.17.0"
        "INFERNAL 1.1.5 (Sep 2023)"                       -> "INFERNAL 1.1.5"
    """
    if not tool_info:
        return ""
    text = str(tool_info).strip().lstrip("#").strip()
    first_line = text.splitlines()[0].strip()
    # "blastx: 2.17.0+" style
    if ":" in first_line and first_line.split(":", 1)[0].strip().lower().startswith(
        "blast"
    ):
        name, _, rest = first_line.partition(":")
        name = (
            name.strip()
            .replace("blastx", "BLASTx")
            .replace("blastn", "BLASTn")
            .replace("blastp", "BLASTp")
        )
        version = rest.strip().split()[0].rstrip("+") if rest.strip() else ""
        return f"{name} {version}".strip()
    # "HMMER 3.4 (...)" / "INFERNAL 1.1.5 (...)" style: take name + first version-ish token
    tokens = first_line.split()
    if not tokens:
        return ""
    name = tokens[0]
    version = tokens[1] if len(tokens) > 1 else ""
    version = version.strip("(),")
    return f"{name} {version}".strip()


# ---------------------------------------------------------------------------
# Data model construction
# ---------------------------------------------------------------------------
def _build_meta(screen: ScreenResult, list_meta: list) -> dict:
    info = getattr(screen, "commec_info", None)
    qinfo = getattr(screen, "query_info", None)
    dbinfo = getattr(screen, "database_info", None)

    tools = ""
    tool_info = getattr(dbinfo, "search_tool_info", None) if dbinfo else None
    if tool_info is not None:
        parts = []
        for attr in (
            "biorisk_search_info",
            "protein_search_info",
            "nucleotide_search_info",
            "low_concern_rna_search_info",
        ):
            stv = getattr(tool_info, attr, None)
            short = _tool_short(getattr(stv, "tool_info", "")) if stv else ""
            if short:
                parts.append(short)
        tools = " · ".join(parts)

    # Database component revisions (helps troubleshoot bug reports).
    revisions = ""
    revs = getattr(dbinfo, "revisions", None) if dbinfo else None
    if isinstance(revs, dict) and revs:
        labels = {
            "biorisk": "biorisk",
            "best_match": "best-match",
            "low_concern": "low-concern",
            "control_lists": "control-lists",
        }
        order = ["biorisk", "best_match", "low_concern", "control_lists"]
        keys = [k for k in order if k in revs] + [k for k in revs if k not in order]
        revisions = " · ".join(labels.get(k, k) + " v" + str(revs[k]) for k in keys)

    file_path = getattr(qinfo, "file", "") if qinfo else ""
    return {
        "file": os.path.basename(file_path) if file_path else "Screening report",
        "nQueries": getattr(qinfo, "number_of_queries", None) if qinfo else None,
        "totalLen": getattr(qinfo, "total_query_length", None) if qinfo else None,
        "version": getattr(info, "commec_version", "") if info else "",
        "output": getattr(info, "json_output_version", "") if info else "",
        "date": getattr(info, "date_run", "") if info else "",
        "time": getattr(info, "time_taken", "") if info else "",
        "tools": tools,
        "revisions": revisions,
    }


def _list_code(includes, acronym: str) -> str:
    """
    Short chip preview for a control list: the two-letter ISO country code when the
    list covers a single country (``includes`` is one alpha-2 code), otherwise the
    list acronym (e.g. groupings like the Australia Group or the EU, which expand to
    many countries).
    """
    if isinstance(includes, str):
        parts = [p.strip() for p in includes.split(",") if p.strip()]
        if len(parts) == 1 and len(parts[0]) == 2 and parts[0].isalpha():
            return parts[0].upper()
    return acronym


def _build_list_meta(screen: ScreenResult) -> tuple:
    """Return (``lists`` for the model, ``name -> {code, acronym, region}`` lookup)."""
    dbinfo = getattr(screen, "database_info", None)
    entries = getattr(dbinfo, "control_list_info", None) if dbinfo else None
    lists = []
    by_name = {}
    for cl in entries or []:
        name = getattr(cl, "name", "") or ""
        acronym = getattr(cl, "acronym", "") or ""
        region = getattr(cl, "region", "") or ""
        status = str(getattr(cl, "status", "") or "")
        code = _list_code(getattr(cl, "includes", ""), acronym)
        lists.append(
            {
                "name": name,
                "code": code,
                "acronym": acronym,
                "region": region,
                "status": status,
            }
        )
        if name:
            by_name[name] = {"code": code, "acronym": acronym, "region": region}
    return lists, by_name


def _best_taxon(taxa: list):
    """Pick the most significant taxon (lowest e-value, tie-break highest identity)."""
    if not taxa:
        return None

    def key(t):
        ev = _clean(t.get("evalue"))
        try:
            ev = float(ev)
        except (TypeError, ValueError):
            ev = float("inf")
        if math.isnan(ev):
            ev = float("inf")
        pid = _clean(t.get("percent_identity")) or 0
        try:
            pid = float(pid)
        except (TypeError, ValueError):
            pid = 0
        return (ev, -pid)

    return sorted(taxa, key=key)[0]


def _target_url(step: str, target: str):
    """
    Build a lookup URL for a best-match target accession: UniProt for protein
    (UniRef) hits, NCBI nuccore for nucleotide hits. Returns None otherwise.
    """
    if not target:
        return None
    if step == ScreenStep.TAXONOMY_AA:
        return "https://www.uniprot.org/uniprotkb/%s/" % target
    if step == ScreenStep.TAXONOMY_NT:
        return "https://www.ncbi.nlm.nih.gov/nuccore/%s/" % target
    return None


def _taxon_lists(taxon: dict, list_lookup: dict) -> list:
    """Join a taxon's ``controlled_by_lists`` (name + source) to run list metadata."""
    out = []
    for entry in taxon.get("controlled_by_lists", []) or []:
        name = entry.get("list", "") or ""
        meta = list_lookup.get(name, {})
        acronym = meta.get("acronym", "")
        out.append(
            {
                "name": name,
                "code": meta.get("code", "") or acronym,
                "acronym": acronym,
                "region": meta.get("region", ""),
                "source": _clean(entry.get("source")) or "",
            }
        )
    return out


def _build_hit(hit, list_lookup: dict) -> dict:
    rec = getattr(hit, "recommendation", None)
    step = str(getattr(rec, "from_step", "") or "")
    region = getattr(hit, "region", None)
    ann = getattr(hit, "annotations", None) or {}

    common = {
        "step": step,
        "rawStatus": str(getattr(rec, "status", "") or ""),
        "name": _clean_text(getattr(hit, "name", "")) or "",
        "desc": _clean_text(getattr(hit, "description", "")) or "",
        "qs": getattr(region, "query_start", None) if region else None,
        "qe": getattr(region, "query_end", None) if region else None,
        "eValue": _fmt_evalue(getattr(region, "e_value", None) if region else None),
        "regulated": None,
        "domain": None,
        "category": None,
        "lists": [],
        "pctId": None,
        "target": None,
        "targetUrl": None,
        "targetDesc": None,
        "taxid": None,
        "genus": None,
        "species": None,
        "nControlled": 0,
        "nNonControlled": 0,
        "coverage": None,
    }

    if step == ScreenStep.BIORISK:
        common["regulated"] = _first(ann.get("regulated"))
        common["domain"] = _first(ann.get("domain"))
        return common

    if step in (
        ScreenStep.LOW_CONCERN_PROTEIN,
        ScreenStep.LOW_CONCERN_RNA,
        ScreenStep.LOW_CONCERN_DNA,
    ):
        common["coverage"] = _clean(ann.get("Coverage: "))
        return common

    if step in (ScreenStep.TAXONOMY_AA, ScreenStep.TAXONOMY_NT):
        common["category"] = _first(ann.get("category"))
        tax = ann.get("controlled_taxonomy") or {}
        stats = tax.get("statistics") or {}
        common["nControlled"] = stats.get("number_of_controlled_taxids", 0) or 0
        common["nNonControlled"] = stats.get("number_of_non-controlled_taxids", 0) or 0

        controlled = tax.get("controlled_taxa") or []
        rep = _best_taxon(controlled) or _best_taxon(
            tax.get("non-controlled_taxa") or []
        )
        if rep:
            common["target"] = _clean(rep.get("target_hit"))
            common["targetUrl"] = _target_url(step, common["target"])
            common["targetDesc"] = _clean_text(rep.get("target_description")) or ""
            common["pctId"] = _fmt_pct(rep.get("percent_identity"))
            common["taxid"] = _clean(rep.get("taxid"))
            common["genus"] = _clean_text(rep.get("genus"))
            common["species"] = _clean_text(rep.get("species"))
            common["lists"] = _taxon_lists(rep, list_lookup)
        return common

    return common


def build_report_model(screen: ScreenResult) -> dict:
    """Flatten a :class:`ScreenResult` into the model consumed by ``report.js``."""
    lists, list_lookup = _build_list_meta(screen)
    meta = _build_meta(screen, lists)

    sequences = []
    for query in getattr(screen, "queries", {}).values():
        status = getattr(query, "status", None)
        sequences.append(
            {
                "id": getattr(query, "query", ""),
                "name": getattr(query, "query", ""),
                "desc": getattr(query, "description", "") or "",
                "length": getattr(query, "length", 0),
                "screenStatus": str(getattr(status, "screen_status", "") or ""),
                "rationale": str(getattr(status, "rationale", "") or ""),
                "hits": [
                    _build_hit(h, list_lookup) for h in getattr(query, "hits", [])
                ],
            }
        )

    return {"meta": meta, "lists": lists, "sequences": sequences}


# ---------------------------------------------------------------------------
# HTML assembly
# ---------------------------------------------------------------------------
def _asset(name: str) -> str:
    return (
        importlib.resources.files("commec")
        .joinpath("utils", "report_assets", name)
        .read_text(encoding="utf-8")
    )


# Open-licensed variable fonts embedded so the report renders identically offline.
# Each is a single variable woff2 covering the whole weight range used by the report.
_FONTS = [
    ("Manrope.woff2", "Manrope", "normal", "400 700"),
    ("CrimsonPro.woff2", "Crimson Pro", "normal", "400 700"),
    ("CrimsonPro-Italic.woff2", "Crimson Pro", "italic", "400 700"),
]


def _font_face_css() -> str:
    """Emit @font-face rules with the woff2 files inlined as base64 data URIs."""
    rules = []
    for filename, family, style, weight in _FONTS:
        data = (
            importlib.resources.files("commec")
            .joinpath("utils", "report_assets", "fonts", filename)
            .read_bytes()
        )
        b64 = base64.b64encode(data).decode("ascii")
        rules.append(
            "@font-face{font-family:'%s';font-style:%s;font-weight:%s;font-display:swap;"
            "src:url(data:font/woff2;base64,%s) format('woff2');}"
            % (family, style, weight, b64)
        )
    return "\n".join(rules)


def _html_escape(text: str) -> str:
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _logo_svg() -> str:
    """Load the masthead logo, strip any XML prolog, and give it a display size."""
    svg = _asset("commec-logo.svg")
    start = svg.find("<svg")
    if start > 0:
        svg = svg[start:]
    opening_tag = svg.split(">", 1)[0]
    if "style=" not in opening_tag:
        svg = svg.replace(
            "<svg", '<svg style="height:24px;width:auto;display:block;"', 1
        )
    return svg.strip()


def render_report_html(screen: ScreenResult) -> str:
    """Build the complete self-contained HTML document for a screen result."""
    model = build_report_model(screen)
    # allow_nan=False is safe because _clean() removed every NaN; escape "</" so a stray
    # "</script>" inside a description can't break out of the data <script> block.
    data_json = json.dumps(model, allow_nan=False, ensure_ascii=False).replace(
        "</", "<\\/"
    )

    css = _font_face_css() + "\n" + _asset("report.css")
    # Inline the logo markup into the renderer (kept as a separate editable asset).
    js = _asset("report.js").replace('"__COMMEC_LOGO__"', json.dumps(_logo_svg()))
    title = _html_escape(
        "Commec Screening Report — " + (model["meta"].get("file") or "")
    )

    return (
        "<!DOCTYPE html>\n"
        '<html lang="en">\n<head>\n'
        '<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        "<title>" + title + "</title>\n"
        "<style>\n" + css + "\n</style>\n"
        "</head>\n<body>\n"
        '<div id="commec-report">'
        '<noscript><div style="max-width:640px;margin:60px auto;padding:0 24px;'
        'font-family:Arial,sans-serif;color:#23285A;line-height:1.6;">'
        "<strong>This Commec screening report needs JavaScript to display.</strong><br>"
        "It is a single self-contained file &mdash; save it and open it in a web browser "
        "(email clients block the embedded script when shown inline).</div></noscript>"
        "</div>\n"
        "<script>window.COMMEC_REPORT = " + data_json + ";</script>\n"
        "<script>\n" + js + "\n</script>\n"
        "</body>\n</html>\n"
    )


def generate_html_from_screen_data(input_data: ScreenResult, output_file: str):
    """
    Render the ScreenResult from Commec Screen as an interactive HTML report.
    Combines all queries from a single run into one HTML document.
    """
    html = render_report_html(input_data)
    output_filename = output_file.strip() + ".html"
    with open(output_filename, "w", encoding="utf-8") as handle:
        handle.write(html)


def generate_html_from_screen_json(input_file: str, output_file: str):
    """Wrapper accepting an input JSON filepath rather than a ScreenResult object."""
    input_data: ScreenResult = get_screen_data_from_json(input_file)
    generate_html_from_screen_data(input_data, output_file)


def main():
    """Convert a JSON output from Commec Screen into an HTML data visualisation."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", dest="in_file", required=True, help="Input json file path"
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="out_file",
        required=True,
        help="Output html filepath, not including .html extension.",
    )
    args = parser.parse_args()
    generate_html_from_screen_json(args.in_file, args.out_file)


if __name__ == "__main__":
    main()
