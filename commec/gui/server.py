"""Tiny web frontend for `commec screen`.

One Flask server, two modes selected by flags:
  - kiosk : bind 127.0.0.1 and auto-open a fullscreen browser
  - LAN   : bind 0.0.0.0 so other machines on the network can connect

Phase 1: FASTA is supplied as a path that exists on the SERVER machine.
Upload and paste-sequence inputs come later.

The server is meant to run inside the `commec-dev` conda env so that the
`commec` binary is on PATH; override with --commec-bin if needed.
"""

import argparse
import atexit
import copy
import csv
import glob
import importlib.resources
import ipaddress
import json
import os
import re
import shlex
import shutil
import signal
import subprocess
import threading
import time
import uuid
import webbrowser
from pathlib import Path

import pycountry
import yaml
from flask import (Flask, Response, jsonify, redirect, request, send_file,
                   session, url_for)
from werkzeug.security import check_password_hash

# CPU/RAM meters read Linux /proc directly (zero-dep; the kiosk path). psutil is
# a purely opportunistic fallback for platforms without /proc (macOS dev/demo):
# never required, never installed by packaging, and never allowed to break the
# server -- if it's absent or misbehaves, the meters just read "--".
try:
    import psutil
except Exception:  # any import failure -> treat as unavailable
    psutil = None

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 200 * 1024 * 1024  # 200 MB upload cap
app.secret_key = os.urandom(32)  # per-process; a restart invalidates sessions
app.config.update(SESSION_COOKIE_HTTPONLY=True, SESSION_COOKIE_SAMESITE="Lax")


@app.errorhandler(413)
def _too_large(_e):
    return jsonify({"error": "Uploaded file is too large (max 200 MB)."}), 413

# ---------------------------------------------------------------------------
# Config, populated from argparse in main(). Kept on the app for test access.
# ---------------------------------------------------------------------------
CFG = {
    "commec_bin": "commec",
    # runs_dir holds one persistent directory per run: the submitted sequence,
    # commec's raw search intermediates, and the polished results all live
    # together and are kept on disk (the box is treated as secured). This
    # enables post-run debugging and future run-resume.
    "runs_dir": Path("runs").resolve(),
    "runs_keep": 0,         # max retained run dirs; 0 = unlimited (forever)
    "default_databases": "",
    "threads": None,        # CPU threads passed to commec (-t); None = don't set
    "browse_root": Path.home(),  # file browser is confined to this directory
    "password_hash": None,  # set from --password-file; None = no auth
    "require_local_auth": False,  # also require the password from localhost
    "presets": [],          # list of {id, label, description, config}
}

# The polished, self-contained artifacts. The GUI's Results panel shows these
# by default; everything else in a run dir (input.fasta, raw search output) is
# tagged "advanced" and revealed only under the Advanced-options toggle.
# {prefix} is filled with the run's output prefix.
PRIMARY_ARTIFACTS = (
    "{prefix}.output.json",
    "{prefix}_summary.html",
    "{prefix}.screen.log",
)


def load_presets(path):
    """Read the presets manifest into a validated list of preset dicts."""
    with open(path, "r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    presets = data.get("presets", [])
    if not isinstance(presets, list) or not presets:
        raise ValueError(f"No presets found in {path}")
    seen = set()
    for p in presets:
        if not all(k in p for k in ("id", "label")):
            raise ValueError(f"Preset missing id/label in {path}: {p!r}")
        if p["id"] in seen:
            raise ValueError(f"Duplicate preset id '{p['id']}' in {path}")
        seen.add(p["id"])
        p.setdefault("description", "")
        p.setdefault("config", {})
        if not isinstance(p["config"], dict):
            raise ValueError(f"Preset '{p['id']}' config must be a mapping")
    return presets


def _preset(preset_id):
    for p in CFG["presets"]:
        if p["id"] == preset_id:
            return p
    return None


def _custom_preset(text):
    """Build a throwaway preset from a raw user-supplied commec config (YAML).
    Returns (preset, None) or (None, error). This is a debug escape hatch: we
    don't validate keys (commec's config merge ignores unknown ones), but we
    block do_cleanup, which crashes commec at end-of-run."""
    if not text.strip():
        return None, "Custom config is empty."
    try:
        cfg = yaml.safe_load(text)
    except yaml.YAMLError as exc:
        return None, f"Custom config isn't valid YAML: {exc}"
    if not isinstance(cfg, dict):
        return None, "Custom config must be YAML key: value pairs (a mapping)."
    if cfg.get("do_cleanup"):
        return None, "Remove do_cleanup -- it crashes commec at end-of-run."
    return {"id": "custom", "label": "custom", "config": cfg}, None


# Top-level database subdirectories under the database dir, with friendly
# labels for messages. Mirrors commec's screen-default-config.yaml layout.
DB_LABELS = {
    "biorisk": "biorisk",
    "low_concern": "low-concern",
    "best_match_protein": "best-match protein (blastx)",
    "best_match_nucleotide": "best-match nucleotide",
    "control_lists": "control lists",
}


def _required_db_dirs(config):
    """Which db subdirectories a preset needs, derived from its skip flags.

    Mirrors commec's gating: best-match protein search runs unless taxonomy is
    skipped (and needs the control lists for regulation calls); best-match
    nucleotide search runs only when neither taxonomy nor nt is skipped.
    """
    dirs = ["biorisk", "low_concern"]
    skip_tx = bool(config.get("skip_taxonomy_search", False))
    skip_nt = bool(config.get("skip_nt_search", False))
    if not skip_tx:
        dirs.append("best_match_protein")
        dirs.append("control_lists")
        if not skip_nt:
            dirs.append("best_match_nucleotide")
    return dirs


def _db_present(path):
    """True if a commec database exists at `path`, which may be a directory
    (control_lists/), an exact file (biorisk/biorisk.hmm), or a BLAST
    file *prefix* (best_match/protein/prot -> prot.00.pin, prot.00.phr, ...)."""
    if not path:
        return False
    return os.path.isdir(path) or os.path.exists(path) or bool(glob.glob(path + "*"))


def _commec_db_paths():
    """Resolve each database to its ACTUAL on-disk path from commec's shipped
    screen-default-config.yaml ({default} expanded). Checking the real configured paths -
    not hardcoded 'best_match' dir names - keeps availability correct whatever DB
    layout the config points at. Falls back to <db_dir>/<name> for any database whose config
    carries no path."""
    db_dir = CFG.get("default_databases") or ""
    base, dbs = db_dir, {}
    try:
        cfg = yaml.safe_load(importlib.resources.files("commec")
                             .joinpath("screen-default-config.yaml").read_text())
        base = (cfg.get("base_paths") or {}).get("default") or db_dir
        dbs = cfg.get("databases", {}) or {}
    except Exception:
        pass
    if base and not base.endswith(os.sep):
        base += os.sep

    def cfg_path(*keys):
        d = dbs
        for k in keys:
            if not isinstance(d, dict):
                return None
            d = d.get(k)
        return d.get("path") if isinstance(d, dict) else None

    def resolve(name, *keys):
        p = cfg_path(*keys)
        return p.replace("{default}", base) if p else os.path.join(db_dir, name)

    return {
        "biorisk":               resolve("biorisk", "biorisk"),
        "low_concern":           resolve("low_concern", "low_concern"),
        "best_match_protein":    resolve("best_match/protein", "best_match", "protein"),
        "best_match_nucleotide": resolve("best_match/nucleotide", "best_match", "nucleotide"),
        "control_lists":         resolve("control_lists", "control_lists"),
    }


def _region_choices(include_single_country_groups=False):
    """Jurisdictions represented by the installed control-list data."""
    control_lists = Path(_commec_db_paths()["control_lists"])
    definitions = control_lists / "region_definitions.json"
    try:
        groups = json.loads(definitions.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        groups = []
    if not isinstance(groups, list):
        groups = []

    choices = []
    group_codes = set()
    memberships = {}
    supported_countries = set()
    for group in groups:
        if not isinstance(group, dict):
            continue
        code = str(group.get("acronym") or "").strip().upper()
        name = str(group.get("name") or "").strip()
        region_codes = list(dict.fromkeys(
            str(region).strip().upper()
            for region in group.get("regions", [])
            if str(region).strip()
        )) if isinstance(group.get("regions"), list) else []
        if not code or not name:
            continue
        group_codes.add(code)
        supported_countries.update(region_codes)
        if region_codes and (len(region_codes) > 1 or include_single_country_groups):
            choices.append({"code": code, "name": name, "group": True})
            for region_code in region_codes:
                memberships.setdefault(region_code, []).append(name)

    for list_info in control_lists.rglob("list_info.csv"):
        try:
            with open(list_info, encoding="utf-8", newline="") as fh:
                supported_countries.update(
                    str(row.get("region_code") or "").strip().upper()
                    for row in csv.DictReader(fh)
                    if str(row.get("region_code") or "").strip()
                )
        except (OSError, csv.Error):
            continue

    countries = {country.alpha_2: country for country in pycountry.countries}
    choices.extend(
        {
            "code": country.alpha_2,
            "name": country.name,
            "group": False,
            "memberships": memberships.get(country.alpha_2, []),
        }
        for country in sorted(
            (
                countries[code]
                for code in supported_countries - group_codes
                if code in countries
            ),
            key=lambda item: item.name,
        )
    )
    return choices


def _missing_databases(preset):
    """Return the required databases that are absent, checked against their real
    configured paths (see _commec_db_paths). [] when no db dir is configured."""
    db_dir = CFG["default_databases"]
    if not db_dir or not os.path.isdir(db_dir):
        return []
    paths = _commec_db_paths()
    return [name for name in _required_db_dirs(preset["config"])
            if not _db_present(paths.get(name, ""))]

# ---------------------------------------------------------------------------
# Offline Plotly. commec's *_summary.html report pulls plotly.js from a CDN
# (include_plotlyjs='cdn'), so it won't render charts on an offline kiosk/LAN.
# We ship a local copy and, when *viewing* a report through this server,
# transparently rewrite the CDN <script> tag to point at our local copy. This
# keeps commec's core untouched; the on-disk HTML still references the CDN.
# ---------------------------------------------------------------------------
VENDOR_DIR = Path(__file__).parent / "vendor"
LOCAL_PLOTLY = "/vendor/plotly.min.js"

# Matches commec's CDN plotly tag at any version. We drop integrity/crossorigin
# along with it: SRI would reject our same-origin copy if bytes ever differ.
_PLOTLY_CDN_RE = re.compile(
    r'<script\b[^>]*\bsrc="https://cdn\.plot\.ly/plotly-[^"]*"[^>]*></script>'
)


def _localize_plotly(html):
    """Point commec's CDN plotly <script> tags at our local copy."""
    return _PLOTLY_CDN_RE.sub(
        f'<script charset="utf-8" src="{LOCAL_PLOTLY}"></script>', html
    )


@app.route("/vendor/plotly.min.js")
def vendor_plotly():
    return send_file(VENDOR_DIR / "plotly-3.0.1.min.js")


ASSETS_DIR = Path(__file__).parent / "assets"


@app.route("/assets/<path:name>")
def asset(name):
    """Serve static brand assets (logo, etc.), confined to assets/."""
    target = (ASSETS_DIR / name).resolve()
    if ASSETS_DIR.resolve() not in target.parents or not target.is_file():
        return jsonify({"error": "not found"}), 404
    return send_file(target)


# ---------------------------------------------------------------------------
# Network addresses. The GUI header shows how to reach this machine. We read
# the live interfaces rather than hard-coding, so a kiosk just works wherever
# it's plugged in. Linux-focused (iproute2 + /sys), used on the Debian kiosk;
# returns [] elsewhere so the header simply hides the section.
# ---------------------------------------------------------------------------
_TAILSCALE_NET = ipaddress.ip_network("100.64.0.0/10")  # Tailscale CGNAT range


def _iface_kind(iface, ip):
    """Classify an interface as ethernet/wifi/tailscale, or None to skip."""
    try:
        if iface.startswith("tailscale") or ipaddress.ip_address(ip) in _TAILSCALE_NET:
            return "tailscale"
    except ValueError:
        pass
    base = f"/sys/class/net/{iface}"
    if os.path.exists(f"{base}/phy80211") or os.path.exists(f"{base}/wireless"):
        return "wifi"
    if os.path.exists(f"{base}/device"):
        return "ethernet"
    return None  # loopback, docker, bridges, etc.


def network_addresses():
    """This host's IPv4 addresses, one per interface, by kind. [] if unknown."""
    try:
        out = subprocess.run(["ip", "-j", "-4", "addr", "show"],
                             capture_output=True, text=True, timeout=3, check=True).stdout
        data = json.loads(out)
    except (OSError, subprocess.SubprocessError, ValueError):
        return []
    rank = {"ethernet": 0, "wifi": 1, "tailscale": 2}
    found = []
    for entry in data:
        iface = entry.get("ifname", "")
        if iface == "lo":
            continue
        for addr in entry.get("addr_info", []):
            ip = addr.get("local")
            if not ip:
                continue
            kind = _iface_kind(iface, ip)
            if kind:
                found.append({"kind": kind, "iface": iface, "ip": ip})
            break  # first IPv4 per interface is enough
    found.sort(key=lambda a: (rank.get(a["kind"], 9), a["iface"]))
    return found


@app.route("/net")
def net():
    return jsonify({"addresses": network_addresses()})


# ---------------------------------------------------------------------------
# System resources for the live meters. Screening is CPU-bound, so the CPU
# meter shows the box working. Linux /proc only; returns null elsewhere.
# ---------------------------------------------------------------------------
_CPU_PREV = {"total": 0, "idle": 0}


def _cpu_percent():
    """System CPU busy % since the previous call (None until the 2nd call)."""
    try:
        vals = [int(x) for x in open("/proc/stat").readline().split()[1:]]
    except (OSError, ValueError):
        try:
            return psutil.cpu_percent(interval=None) if psutil else None
        except Exception:
            return None
    idle = vals[3] + (vals[4] if len(vals) > 4 else 0)  # idle + iowait
    total = sum(vals)
    prev_total, prev_idle = _CPU_PREV["total"], _CPU_PREV["idle"]
    _CPU_PREV["total"], _CPU_PREV["idle"] = total, idle
    dt = total - prev_total
    if prev_total == 0 or dt <= 0:
        return None  # need a baseline / no elapsed ticks yet
    return max(0.0, min(100.0, 100.0 * (1 - (idle - prev_idle) / dt)))


def _mem_info():
    """(percent_used, used_bytes, total_bytes) from /proc/meminfo, or None."""
    try:
        info = {}
        for line in open("/proc/meminfo"):
            key, _, val = line.partition(":")
            info[key] = int(val.split()[0]) * 1024  # kB -> bytes
    except (OSError, ValueError):
        try:
            if psutil:
                m = psutil.virtual_memory()
                return (m.percent, m.used, m.total)
        except Exception:
            pass
        return None
    total = info.get("MemTotal")
    avail = info.get("MemAvailable", info.get("MemFree", 0))
    if not total:
        return None
    used = total - avail
    return (100.0 * used / total, used, total)


def _disk_info():
    """Disk usage of the filesystem holding the runs dir, or None. This is the
    disk that fills up as runs accumulate, so it's the one worth showing."""
    try:
        u = shutil.disk_usage(CFG["runs_dir"])
    except OSError:
        return None
    pct = (100.0 * u.used / u.total) if u.total else None
    return {"pct": pct, "free": u.free, "used": u.used, "total": u.total}


def _dir_size(d):
    """Total bytes of all files under a run dir (intermediates included)."""
    total = 0
    for p in d.rglob("*"):
        try:
            if p.is_file():
                total += p.stat().st_size
        except OSError:
            pass
    return total


@app.route("/metrics")
def metrics():
    cpu = _cpu_percent()
    mem = _mem_info()
    return jsonify({
        "cpu": cpu,
        "ram": ({"pct": mem[0], "used": mem[1], "total": mem[2]} if mem else None),
        "disk": _disk_info(),
    })


# ---------------------------------------------------------------------------
# Server-side file browser for the "Server file path" tab. Confined to
# CFG["browse_root"] so clients can only see/pick files under it (symlinks that
# escape resolve outside and get clamped back to the root).
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Removable media (USB). A background thread mounts removable partitions
# read-only via udisksctl (no root; see deploy-notes.md for the host grant);
# the page surfaces mounted volumes as tabs. Read-only access only.
# ---------------------------------------------------------------------------
def _lsblk():
    """Flat list of block-device nodes from lsblk -J, or []."""
    try:
        out = subprocess.run(
            ["lsblk", "-J", "-o", "NAME,PATH,RM,TYPE,FSTYPE,MOUNTPOINT,LABEL"],
            capture_output=True, text=True, timeout=5, check=True).stdout
        data = json.loads(out)
    except (OSError, subprocess.SubprocessError, ValueError):
        return []
    flat = []

    def walk(nodes):
        for n in nodes:
            flat.append(n)
            walk(n.get("children", []))

    walk(data.get("blockdevices", []))
    return flat


def _rm(node):
    return str(node.get("rm")).lower() in ("1", "true")


def _mountpoint(node):
    mp = node.get("mountpoint")
    if mp:
        return mp
    return next((m for m in (node.get("mountpoints") or []) if m), None)


def _removable_volumes():
    """Mounted removable partitions: [{label, mountpoint, dev}]."""
    vols = []
    for n in _lsblk():
        if n.get("type") == "part" and _rm(n) and n.get("fstype"):
            mp = _mountpoint(n)
            if mp:
                vols.append({"label": n.get("label") or n.get("name"),
                             "mountpoint": mp, "dev": n.get("path")})
    return vols


def _automount_removable():
    """Mount any unmounted removable partition read-only (best-effort)."""
    for n in _lsblk():
        if (n.get("type") == "part" and _rm(n) and n.get("fstype")
                and not _mountpoint(n) and n.get("path")):
            try:
                subprocess.run(["udisksctl", "mount", "-b", n["path"], "-o", "ro"],
                               capture_output=True, text=True, timeout=20)
            except (OSError, subprocess.SubprocessError):
                pass


def _automount_loop(interval):
    while True:
        try:
            _automount_removable()
        except OSError:
            pass
        time.sleep(interval)


@app.route("/media")
def media():
    return jsonify({"volumes": _removable_volumes()})


def _allowed_roots():
    """Directories the browser may show: the configured root + media mounts."""
    roots = [Path(CFG["browse_root"]).resolve()]
    for v in _removable_volumes():
        try:
            roots.append(Path(v["mountpoint"]).resolve())
        except (OSError, ValueError):
            pass
    return roots


@app.route("/browse")
def browse():
    roots = _allowed_roots()
    try:
        cur = Path(request.args.get("path") or roots[0]).resolve()
    except (OSError, ValueError):
        cur = roots[0]
    # Confine to whichever allowed root contains cur (USB mount or home).
    base = next((r for r in roots if cur == r or r in cur.parents), None)
    if base is None or not cur.is_dir():
        cur = base = roots[0]
    entries = []
    try:
        for p in sorted(cur.iterdir(), key=lambda x: (not x.is_dir(), x.name.lower())):
            if p.name.startswith("."):
                continue
            try:
                is_dir = p.is_dir()
            except OSError:
                continue
            entries.append({
                "name": p.name, "path": str(p), "dir": is_dir,
                "size": (p.stat().st_size if not is_dir else None),
            })
    except OSError as exc:
        return jsonify({"error": f"Cannot read directory: {exc}"}), 400
    return jsonify({
        "path": str(cur),
        "parent": (str(cur.parent) if cur != base else None),
        "entries": entries,
    })


# ---------------------------------------------------------------------------
# Job tracking. One job at a time behind a lock: screening is CPU-bound and
# --threads already saturates the box, so a queue would be over-engineering.
# ---------------------------------------------------------------------------
JOBS = {}          # job_id -> job dict
JOB_LOCK = threading.Lock()


# ---------------------------------------------------------------------------
# Input handling. The browser can submit a server-side FASTA path or pasted
# sequence text (FASTA, or rows copied from a spreadsheet). Both are normalised
# into a FASTA file in the run dir.
# ---------------------------------------------------------------------------
def _clean_seq(s):
    # Detection heuristic ONLY (see _looks_like_dna). Pasted/uploaded sequences are NOT sanitised
    # before screening: raw content goes to commec as-is, so a malformed sequence fails loudly at
    # commec rather than being silently repaired here.
    return re.sub(r"[^A-Za-z]", "", s)


def _looks_like_dna(s):
    s = _clean_seq(s)
    if len(s) < 10:
        return False
    acgtn = sum(c in "ACGTNUacgtnu" for c in s)
    return acgtn / len(s) >= 0.9


def _parse_fasta_text(text):
    """Parse FASTA text into [(name, seq)] records."""
    records, name, seq = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq)))
            name, seq = line[1:].strip(), []
        elif line.strip():
            seq.append(line.strip())
    if name is not None:
        records.append((name, "".join(seq)))
    return records


def _parse_spreadsheet_text(text):
    """Tab/comma-delimited rows (Excel/CSV paste): the most DNA-like field is
    the sequence; another field, if present, names it. Rows with no DNA-like
    field (e.g. a header row) are skipped."""
    records = []
    for i, line in enumerate(text.splitlines()):
        fields = [f.strip() for f in re.split(r"[\t,]", line) if f.strip()]
        dna = [f for f in fields if _looks_like_dna(f)]
        if not dna:
            continue
        seq = max(dna, key=len)
        names = [f for f in fields if f is not seq]
        records.append((names[0] if names else f"sequence_{i + 1}", seq))
    return records


def _parse_blocks_text(text):
    """Header-less raw sequences: blocks separated by a blank line, each block's
    lines joined into one sequence (so a single sequence may wrap over lines)."""
    records = []
    for i, block in enumerate(re.split(r"\n\s*\n", text.strip())):
        # Join the block's (possibly wrapped) lines into one sequence WITHOUT sanitising content.
        seq = "".join(line.strip() for line in block.splitlines())
        if seq:
            records.append((f"sequence_{i + 1}", seq))
    return records


def parse_pasted_sequences(text):
    """Turn pasted text into FASTA records, accepting FASTA, raw sequences
    (separated by blank lines), or spreadsheet rows."""
    text = text.strip()
    if not text:
        return []
    if re.search(r"^\s*>", text, re.M):
        return _parse_fasta_text(text)
    if "\t" in text or "," in text:
        return _parse_spreadsheet_text(text)
    return _parse_blocks_text(text)


def _normalise_records(records):
    """Drop empty sequences; make names FASTA-safe and unique."""
    out, used = [], set()
    for i, (name, seq) in enumerate(records):
        if not seq:
            continue
        name = re.sub(r"\s+", "_", (name or "").strip())
        name = re.sub(r"[<>&\"'`]", "", name) or f"sequence_{i + 1}"  # keep names HTML-safe
        base, n = name, 2
        while name in used:
            name, n = f"{base}_{n}", n + 1
        used.add(name)
        out.append((name, seq))
    return out


def _write_fasta(records, path):
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for j in range(0, len(seq), 70):
                fh.write(seq[j:j + 70] + "\n")


def _build_command(fasta, outdir, config_path, prefix, resume=False):
    """Build the commec screen argv list for a run.

    Args are always passed as a list (never shell=True) so user-supplied
    paths can never be interpreted as shell syntax. The chosen config is
    passed via -y; the database directory is injected via -d so it overrides
    whatever (if anything) the config declares.

    resume=True adds -R, telling commec to re-use pre-existing output in the
    dir (used to continue an interrupted run -- see /runs/<id>/resume).
    """
    cmd = [CFG["commec_bin"], "screen", "-y", str(config_path)]
    if CFG["default_databases"]:
        cmd += ["-d", CFG["default_databases"]]
    if CFG["threads"]:
        cmd += ["-t", str(CFG["threads"])]
    if resume:
        cmd += ["-R"]
    # -o as a basename (outdir/prefix, no trailing sep) names the outputs after
    # the run label, independent of the input filename.
    cmd += ["-o", str(outdir / prefix), fasta]
    return cmd


# User-fixable commec errors we lift out of the (advanced-only) terminal log
# into a plain message on the results card. These are caused by the operator's
# input -- unlike internal failures -- so they deserve to be seen directly.
def _error_hint(log_lines):
    """A plain, actionable message for a user-caused error in commec's output,
    or None. Add cases here as more get reported from real usage."""
    text = "\n".join(log_lines)
    if "Duplicate sequence identifier generated" in text:
        m = re.search(r'Duplicate sequence identifier generated:\s*"([^"]+)"', text)
        which = f' (e.g. "{m.group(1)}")' if m else ""
        return ("Duplicate sequence names in your input" + which + ". Each FASTA "
                "header must be unique within its first 64 characters -- rename "
                "the duplicates and run again.")
    return None


def _job_meta(job, status, summary=None):
    """The meta.json payload for a job at a given status (the durable marker
    that survives a restart and drives the results/recent views)."""
    return {
        "id": job["id"], "label": job["label"], "prefix": job["prefix"],
        "fasta": job.get("fasta"), "status": status,
        "returncode": job["returncode"], "created": job["created"],
        "finished": job.get("finished"), "summary": summary,
        "error_hint": job.get("error_hint"),
    }


def _save_meta(rundir, meta):
    try:
        (rundir / "meta.json").write_text(json.dumps(meta), encoding="utf-8")
        return True
    except OSError:
        return False


def _finalize_job(job):
    """Mark a run finished in place. The run dir (sequence + intermediates +
    results) is kept; nothing is copied or wiped. Overwrites the start marker
    with the terminal status + parsed verdict."""
    rundir = job["dir"]
    job["finished"] = time.time()
    meta = _job_meta(job, job["status"], summary=_summarize_output(rundir))
    if not _save_meta(rundir, meta):
        job["log"].append("[gui] WARNING: could not write run marker.")
    _sweep_runs()  # apply any rolling-retention cap


def _run_job(job):
    """Spawn commec and pump its merged stdout/stderr into the job log."""
    cmd = job["cmd"]
    try:
        # start_new_session: commec leads its own process group, so we can kill
        # the whole tree (commec + blast/hmmer children) on cancel/shutdown/sweep.
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            start_new_session=True,
        )
    except FileNotFoundError:
        job["log"].append(f"ERROR: commec binary not found: {cmd[0]}")
        job["status"] = "error"
        job["returncode"] = 127
        _finalize_job(job)
        job["done"] = True
        return

    job["proc"] = proc
    job["pgid"] = proc.pid  # session leader: pgid == pid
    try:
        (job["dir"] / ".pgid").write_text(str(proc.pid), encoding="utf-8")
    except OSError:
        pass
    job["status"] = "running"
    # Durable "running" marker: if the process is later killed (crash/restart)
    # this is what the sweep promotes to "interrupted", flagging it resumable.
    _save_meta(job["dir"], _job_meta(job, "running"))
    for line in proc.stdout:
        job["log"].append(line.rstrip("\n"))
    proc.wait()
    job["returncode"] = proc.returncode
    if job.get("stopped"):
        # Deliberately halted via /stop: treat like a crash -- resumable, not a
        # dead-end error (the partial step is dropped on resume).
        job["status"] = "interrupted"
    else:
        job["status"] = "done" if proc.returncode == 0 else "error"
    if job["status"] == "error":
        job["error_hint"] = _error_hint(job["log"])
    # Write the completion marker BEFORE marking done, so a client that loads
    # results on the 'end' event sees a finalized run.
    _finalize_job(job)
    job["done"] = True


# ---------------------------------------------------------------------------
# Cleanup: kill process groups orphaned by crashes/restarts, and apply the
# optional rolling-retention cap. Run dirs are otherwise kept (never wiped).
# ---------------------------------------------------------------------------
def _kill_pgid(pgid):
    try:
        os.killpg(int(pgid), signal.SIGKILL)
    except (ProcessLookupError, ValueError, OSError):
        pass


_TERMINAL_STATUS = ("done", "error", "interrupted")


def _sweep_orphans():
    """Promote runs abandoned by a server crash/restart to 'interrupted' (which
    flags them resumable), killing any process group still orphaned from the
    dead run. Never deletes anything; runs with a terminal status (cleanly
    finished, errored, or already flagged) are left untouched.

    Note: a run killed while the server is still alive finalizes itself as
    'error' through _run_job -- 'interrupted' only arises when the server
    process itself died mid-run, so its marker is stuck at 'running'."""
    runs = CFG["runs_dir"]
    if not runs.is_dir():
        return
    with JOB_LOCK:
        live = {j["dir"].name for j in JOBS.values() if not j["done"]}
    for d in runs.iterdir():
        if not d.is_dir() or d.name in live:
            continue
        meta = _read_meta(d)
        if meta.get("status") in _TERMINAL_STATUS:
            continue
        # Abandoned mid-run: reap any lingering process group, then flag it. We
        # unlink .pgid so a recycled PID can't be re-killed on a later sweep.
        pgid_file = d / ".pgid"
        if pgid_file.is_file():
            try:
                _kill_pgid(pgid_file.read_text(encoding="utf-8").strip())
                pgid_file.unlink()
            except OSError:
                pass
        meta.update(id=meta.get("id", d.name), label=meta.get("label", d.name),
                    status="interrupted")
        meta.setdefault("finished", d.stat().st_mtime)
        _save_meta(d, meta)


def _drop_interrupted_artifact(d, prefix):
    """Delete the newest file under input_<prefix>/ or output_<prefix>/ -- the
    step commec was writing when it was killed. Steps run sequentially, so the
    newest-mtime artifact is the interrupted (possibly truncated) one; removing
    it forces just that step to re-run while -R reuses the older, complete ones.
    (.screen.log is excluded by only scanning the two step subdirs.) Returns the
    deleted file's path relative to the run dir, or None if there's nothing."""
    cands = []
    for sub in (d / f"input_{prefix}", d / f"output_{prefix}"):
        if sub.is_dir():
            cands += [p for p in sub.rglob("*") if p.is_file()]
    if not cands:
        return None
    newest = max(cands, key=lambda p: p.stat().st_mtime)
    try:
        rel = newest.relative_to(d).as_posix()
        newest.unlink()
        return rel
    except OSError:
        return None


def _sweep_runs():
    """Apply the rolling-retention cap, if any (default: keep everything).
    Deletes whole run dirs (intermediates included) oldest-first past the cap."""
    keep = CFG["runs_keep"]
    if not keep or keep <= 0 or not CFG["runs_dir"].is_dir():
        return
    dirs = [d for d in CFG["runs_dir"].iterdir() if d.is_dir()]
    dirs.sort(key=lambda d: d.stat().st_mtime, reverse=True)
    for d in dirs[keep:]:
        shutil.rmtree(d, ignore_errors=True)


def _sweep():
    _sweep_orphans()
    _sweep_runs()


def _sweep_loop(interval):
    while True:
        time.sleep(interval)
        try:
            _sweep()
        except OSError:
            pass


def _kill_running_jobs():
    """Kill the process group of any in-flight run (server shutdown)."""
    for job in list(JOBS.values()):
        if not job["done"] and job.get("pgid"):
            try:
                os.killpg(job["pgid"], signal.SIGTERM)
            except OSError:
                pass


# ---------------------------------------------------------------------------
# Authentication. The walk-up kiosk (127.0.0.1) may run without a hash because
# physical presence is the trust. A non-localhost client is never allowed in
# until a password hash exists: otherwise a failed first-run password setup
# would silently expose the LAN GUI. Once a hash exists, every non-localhost
# request needs a valid session (and --require-local-auth includes the kiosk).
# Only meaningful over HTTPS (login/cookie would otherwise be cleartext).
# ---------------------------------------------------------------------------
_LOCAL_ADDRS = {"127.0.0.1", "::1", "localhost"}
# Endpoints reachable without a session (login itself + brand assets it shows).
_AUTH_EXEMPT = {"login", "logout", "asset"}


def _auth_required():
    if not CFG["password_hash"]:
        return False
    if not CFG["require_local_auth"] and request.remote_addr in _LOCAL_ADDRS:
        return False
    return True


def _safe_next(target):
    """Only allow same-site relative redirects (no open-redirect)."""
    if target and target.startswith("/") and not target.startswith("//"):
        return target
    return "/"


@app.before_request
def _gate():
    # Do this before the login/assets exemption: there is no credential a remote
    # client could submit until first-run setup has written a hash, so expose no
    # part of the GUI over LAN in that state.
    if not CFG["password_hash"] and request.remote_addr not in _LOCAL_ADDRS:
        return jsonify({"error": "Remote access is unavailable until an operator password is set."}), 503
    if request.endpoint in _AUTH_EXEMPT:
        return
    if _auth_required() and not session.get("authed"):
        if request.method == "GET" and "text/html" in request.headers.get("Accept", ""):
            return redirect(url_for("login", next=request.full_path.rstrip("?")))
        return jsonify({"error": "Authentication required."}), 401


@app.route("/login", methods=["GET", "POST"])
def login():
    if request.method == "POST":
        password = request.form.get("password") or ""
        nxt = _safe_next(request.form.get("next"))
        if CFG["password_hash"] and check_password_hash(CFG["password_hash"], password):
            session["authed"] = True
            return redirect(nxt)
        return redirect(url_for("login", next=nxt, error=1))
    # Already authed (or auth off) -> no reason to show the login screen.
    if not _auth_required() or session.get("authed"):
        return redirect(_safe_next(request.args.get("next")))
    resp = send_file(Path(__file__).parent / "login.html")
    resp.headers["Cache-Control"] = "no-store"
    return resp


@app.route("/logout")
def logout():
    session.clear()
    return redirect(url_for("login"))


@app.route("/")
def index():
    # no-store so kiosks/clients always fetch the current page (no stale UI).
    resp = send_file(Path(__file__).parent / "index.html")
    resp.headers["Cache-Control"] = "no-store"
    return resp


@app.route("/config")
def config():
    """Expose the preset and regional choices the page renders."""
    presets = []
    for p in CFG["presets"]:
        missing = _missing_databases(p)
        presets.append({
            "id": p["id"],
            "label": p["label"],
            "config": p["config"],
            "recommended": p.get("recommended"),
            "available": not missing,
            "missing": [DB_LABELS.get(m, m) for m in missing],
        })
    return jsonify({"presets": presets, "regions": _region_choices()})


def _label_exists(label):
    """True if a retained run already uses this label (keeps labels unique)."""
    rdir = CFG["runs_dir"]
    if not rdir.is_dir():
        return False
    for d in rdir.iterdir():
        if d.is_dir() and _read_meta(d).get("label") == label:
            return True
    return False


@app.route("/screen", methods=["POST"])
def screen():
    preset_id = (request.form.get("preset") or "").strip()
    if preset_id == "__custom__":
        preset, err = _custom_preset(request.form.get("custom_config") or "")
        if err:
            return jsonify({"error": err}), 400
    else:
        preset = _preset(preset_id)
        if preset is None:
            return jsonify({"error": f"Unknown preset: {preset_id!r}"}), 400

    missing = _missing_databases(preset)
    if missing:
        labels = ", ".join(DB_LABELS.get(m, m) for m in missing)
        return jsonify({"error":
            f"The '{preset['label']}' preset needs databases not installed on "
            f"this server: {labels}. Install them with `commec setup`, or pick "
            f"a preset that doesn't require them."}), 400

    raw_regions = (request.form.get("regions") or "all").strip()
    if raw_regions.lower() == "all":
        regions = "all"
    else:
        region_codes = list(dict.fromkeys(
            code.strip().upper() for code in raw_regions.split(",") if code.strip()
        ))
        allowed_regions = {
            choice["code"]
            for choice in _region_choices(include_single_country_groups=True)
        }
        unknown_regions = [code for code in region_codes if code not in allowed_regions]
        if not region_codes or unknown_regions:
            detail = ", ".join(unknown_regions) or raw_regions
            return jsonify({"error": f"Unknown regulatory jurisdiction: {detail}."}), 400
        regions = ",".join(region_codes)

    run_config = copy.deepcopy(preset["config"])
    databases_config = run_config.setdefault("databases", {})
    if not isinstance(databases_config, dict):
        return jsonify({"error": "Preset databases config must be a mapping."}), 400
    control_list_config = databases_config.setdefault("control_lists", {})
    if not isinstance(control_list_config, dict):
        return jsonify({"error": "Preset control-list config must be a mapping."}), 400
    control_list_config["regions"] = regions

    # Required run label -> output-file prefix, sanitised to a filename-safe
    # token (timestamp fallback only if the label is all punctuation).
    label = (request.form.get("label") or "").strip()
    if not label:
        return jsonify({"error": "A run label is required."}), 400
    if _label_exists(label):
        return jsonify({"error":
            "Can't accept run label: run label already chosen."}), 400
    prefix = re.sub(r"[^A-Za-z0-9._-]+", "_", label).strip("._-")
    if not prefix:
        prefix = time.strftime("run-%Y%m%d-%H%M%S")

    # Screen EVERY provided input source as one batch (paste + upload + a picked
    # server file), rather than a single "mode" -- no queue; the user merges
    # inputs by filling in more than one. Silently screening a few extra
    # sequences is a safer failure than missing some, so we combine. All input
    # is materialised into one input.fasta INSIDE the run dir, which is required,
    # not just tidy: commec writes its output next to the INPUT file (it ignores
    # the directory part of -o; see screen_io._get_output_prefixes), so the input
    # must live in the run dir or the output scatters into the source directory.
    records, sources = [], []  # sources: per-source "(N from X)" for the run note

    text = request.form.get("sequence_text") or ""
    if text.strip():
        r = parse_pasted_sequences(text)
        if not r:
            return jsonify({"error": "Couldn't find any sequences in the pasted "
                "text. Paste FASTA, or spreadsheet rows with a sequence column."}), 400
        records += r
        sources.append(f"{len(r)} pasted")

    f = request.files.get("upload_file")
    if f and f.filename:
        r = parse_pasted_sequences(f.read().decode("utf-8", errors="replace"))
        if not r:
            return jsonify({"error": f"No sequences found in the uploaded file "
                f"'{f.filename}' (expected FASTA, or CSV/TSV with a sequence "
                "column)."}), 400
        records += r
        sources.append(f"{len(r)} from {f.filename}")

    picked = (request.form.get("fasta_file") or "").strip()
    if picked:
        # Confine to the same roots the browser exposes (home + USB mounts), so
        # a hand-crafted POST can't screen arbitrary server files. Resolve first
        # to defeat symlink/.. escapes, matching /browse and /results/file.
        target = Path(picked).resolve()
        if not any(target == r or r in target.parents for r in _allowed_roots()):
            return jsonify({"error": "That file is outside the allowed directory."}), 400
        if not target.is_file():
            return jsonify({"error": f"FASTA path not found on server: {target}"}), 400
        try:
            r = _parse_fasta_text(target.read_text(encoding="utf-8"))
        except OSError as exc:
            return jsonify({"error": f"Could not read {target}: {exc}"}), 400
        if not r:
            return jsonify({"error": f"No FASTA records found in {target.name}."}), 400
        records += r
        sources.append(f"{len(r)} from {target.name}")

    if not records:
        return jsonify({"error": "No sequences to screen. Paste, upload, or pick "
                                 "a file (you can combine them)."}), 400

    note = []
    if len(sources) > 1:
        note.append("Combined inputs: " + ", ".join(sources) + ".")
    note.append(
        "Regulatory jurisdiction: "
        + ("all jurisdictions." if regions == "all" else regions + ".")
    )
    # Normalise the COMBINED list so names are made unique across all sources.
    records = _normalise_records(records)

    with JOB_LOCK:
        busy = [j for j in JOBS.values() if not j["done"]]
        if busy:
            return jsonify({"error": "A screen is already running. "
                                     "Wait for it to finish."}), 409

        job_id = uuid.uuid4().hex[:12]
        rundir = CFG["runs_dir"] / job_id  # persistent: sequence + intermediates + results
        rundir.mkdir(parents=True, exist_ok=True)

        # Write the preset's config into the run dir (also a reproducibility
        # record, kept with the results) and point commec at it.
        config_path = rundir / "config.used.yaml"
        with open(config_path, "w", encoding="utf-8") as fh:
            yaml.safe_dump(run_config, fh, default_flow_style=False)

        fasta = str(rundir / "input.fasta")
        _write_fasta(records, fasta)
        note.insert(0, f"Prepared {len(records)} sequence(s) for screening.")

        cmd = _build_command(fasta, rundir, config_path, prefix)
        job = {
            "id": job_id,
            "label": label or prefix,
            "prefix": prefix,
            "fasta": fasta,
            "created": time.time(),
            "finished": None,
            "cmd": cmd,
            "dir": rundir,
            "log": [],
            "status": "starting",
            "returncode": None,
            "done": False,
            "stopped": False,
            "proc": None,
            "pgid": None,
        }
        JOBS[job_id] = job

    for line in note:
        job["log"].append("[gui] " + line)
    job["log"].append("$ " + " ".join(shlex.quote(c) for c in cmd))
    threading.Thread(target=_run_job, args=(job,), daemon=True).start()
    return jsonify({"job_id": job_id})


@app.route("/status")
def status():
    """Current screener state, shared across clients so a second browser can
    pick up an in-progress run. 'active' is the running job, or the most recent
    one if idle, so a late-joining client can attach to its log and results."""
    with JOB_LOCK:
        running = [j for j in JOBS.values() if not j["done"]]
        if running:
            job, busy = running[0], True
        else:
            job = max(JOBS.values(), key=lambda j: j["created"], default=None)
            busy = False
    active = None if job is None else {
        "job_id": job["id"],
        "label": job["label"],
        "status": job["status"],
        "returncode": job["returncode"],
    }
    return jsonify({"busy": busy, "active": active})


@app.route("/events/<job_id>")
def events(job_id):
    """Stream the job log live via Server-Sent Events, replaying from start
    so late joins and reconnects see the full history."""
    job = JOBS.get(job_id)
    if not job:
        return jsonify({"error": "unknown job"}), 404

    def stream():
        sent = 0
        while True:
            log = job["log"]
            while sent < len(log):
                payload = json.dumps({"line": log[sent]})
                yield f"event: log\ndata: {payload}\n\n"
                sent += 1
            if job["done"] and sent >= len(job["log"]):
                payload = json.dumps({
                    "status": job["status"],
                    "returncode": job["returncode"],
                    "error_hint": job.get("error_hint"),
                })
                yield f"event: end\ndata: {payload}\n\n"
                return
            time.sleep(0.3)

    return Response(stream(), mimetype="text/event-stream")


def _summarize_output(outdir):
    """Per-query verdicts read from commec's *.output.json, or None.

    Parses the JSON directly (commec isn't importable from the GUI env), so we
    can show a verdict table without re-running `commec flag`.
    """
    jsons = sorted(outdir.glob("*.output.json"))
    if not jsons:
        return None
    try:
        data = json.loads(jsons[0].read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    # Overall verdict = the most severe per-query status. Skip ranks just below Flag: an
    # unscreened query (e.g. too long) is treated as nearly as serious as a flagged one, and
    # above Warning/Pass. Matched by base word so commec's variants ("Skip (too long)",
    # "Pass (Skipped Taxonomy)") are counted; unknown/error defaults to the Skip level.
    sev = {"Flag": 4, "Skip": 3, "Warning": 2, "Pass": 1}
    def _sev(s):
        s = str(s or "")
        return next((n for base, n in sev.items() if s.startswith(base)), 3)
    rows = []
    for name, v in (data.get("queries") or {}).items():
        st = v.get("status") or {}
        rows.append({
            "name": name,
            "length": v.get("length"),
            "status": st.get("screen_status", "Unknown"),
            "rationale": st.get("rationale", "") or "",
        })
    if not rows:
        return None
    overall = max((r["status"] for r in rows), key=_sev)
    return {"n": len(rows), "overall": overall, "queries": rows}


def _job_dir(job_id):
    """The single persistent directory holding a run's files (same location
    live or finished). Returns None if it doesn't exist."""
    d = CFG["runs_dir"] / job_id
    return d if d.is_dir() else None


def _is_primary(rel):
    """True if a run-dir-relative path is one of the polished PRIMARY_ARTIFACTS
    (matched at the top level, prefix-agnostic via the template suffix)."""
    if rel.parent != Path("."):
        return False  # anything in a subdir (output_*/...) is an intermediate
    for tmpl in PRIMARY_ARTIFACTS:
        if "{prefix}" in tmpl:
            if rel.name.endswith(tmpl.replace("{prefix}", "")):
                return True
        elif rel.name == tmpl:
            return True
    return False


def _list_run_files(d):
    """All files under a run dir, recursively (intermediates included), tagged
    primary vs advanced. Skips dotfiles and the meta.json marker. Names are
    paths relative to the run dir (the file route accepts nested paths)."""
    files = []
    for p in sorted(d.rglob("*")):
        if not p.is_file() or p.name == "meta.json":
            continue
        rel = p.relative_to(d)
        if any(part.startswith(".") for part in rel.parts):
            continue  # dotfiles (.pgid) and anything under a dotdir
        files.append({
            "name": rel.as_posix(),
            "size": p.stat().st_size,
            "advanced": not _is_primary(rel),
        })
    files.sort(key=lambda f: (f["advanced"], f["name"]))
    return files


@app.route("/results/<job_id>")
def results(job_id):
    job = JOBS.get(job_id)
    d = _job_dir(job_id)
    if d is None and not job:
        return jsonify({"error": "unknown job"}), 404
    files = _list_run_files(d) if d is not None else []
    if job:
        status, returncode, label = job["status"], job["returncode"], job["label"]
        error_hint = job.get("error_hint")
    else:
        meta = _read_meta(d)
        status = _effective_status(meta, False)
        returncode = meta.get("returncode")
        label = meta.get("label")
        error_hint = meta.get("error_hint")
    # Resumable whenever the run is interrupted and still has its config -- same
    # rule as /runs. (A live or finished job has status running/done, so this
    # naturally excludes them; being in JOBS is irrelevant -- a stopped run
    # stays in JOBS as 'interrupted' and is very much resumable.)
    resumable = (status == "interrupted"
                 and d is not None and (d / "config.used.yaml").is_file())
    return jsonify({
        "status": status,
        "returncode": returncode,
        "label": label,
        "files": files,
        "summary": _summarize_output(d) if d is not None else None,
        "resumable": resumable,
        "error_hint": error_hint,
    })


def _read_meta(d):
    try:
        return json.loads((d / "meta.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def _effective_status(meta, is_live):
    """Display status for a run on disk. A non-live run whose marker is still
    non-terminal (no marker, 'starting', or 'running') was abandoned by a
    crash/restart -> report it 'interrupted' immediately, without waiting for
    the periodic sweep to rewrite the marker."""
    s = meta.get("status")
    if is_live:
        return s or "running"
    return s if s in _TERMINAL_STATUS else "interrupted"


@app.route("/results/<job_id>/file/<path:name>")
def result_file(job_id, name):
    base = _job_dir(job_id)
    if base is None:
        return jsonify({"error": "not found"}), 404
    # Resolve and confine to the run's directory (no path traversal).
    target = (base / name).resolve()
    if base.resolve() not in target.parents or not target.is_file():
        return jsonify({"error": "not found"}), 404
    view = request.args.get("view") == "1"
    # When viewing an HTML report in-browser, swap commec's CDN plotly tag for
    # our local copy so charts render offline. Downloads keep the file as-is.
    if view and target.suffix.lower() in (".html", ".htm"):
        html = _localize_plotly(target.read_text(encoding="utf-8"))
        return Response(html, mimetype="text/html")
    return send_file(target, as_attachment=not view)


@app.route("/runs")
def runs():
    """Finished runs (newest first) for the Recent results panel, read from the
    persisted meta.json markers so the list survives restarts. Currently-live
    runs are shown by the live Run panel, so they're skipped here."""
    rdir = CFG["runs_dir"]
    with JOB_LOCK:
        live = {j["dir"].name for j in JOBS.values() if not j["done"]}
    entries = []
    if rdir.is_dir():
        for d in rdir.iterdir():
            if not d.is_dir() or d.name in live:
                continue
            meta = _read_meta(d)
            status = _effective_status(meta, False)
            finished = meta.get("finished") or d.stat().st_mtime
            summary = meta.get("summary") or {}
            # The verdict is cached in the marker at finalize, but a run may have
            # finished without it (e.g. a resume whose marker predates the new
            # output.json). Fall back to deriving it live, and self-heal the
            # marker so the list and the detail view never disagree.
            if not summary and status in ("done", "error"):
                summary = _summarize_output(d) or {}
                if summary:
                    meta["summary"] = summary
                    _save_meta(d, meta)
            entries.append((finished, {
                "id": meta.get("id", d.name),
                "label": meta.get("label", d.name),
                "status": status,
                "returncode": meta.get("returncode"),
                "finished": finished,
                "overall": summary.get("overall"),
                "n": summary.get("n"),
                "size": _dir_size(d),  # bytes on disk (intermediates included)
                # Resumable only if we still have what -R needs (the config; the
                # input is re-checked at resume time, since it may be external).
                "resumable": (status == "interrupted"
                              and (d / "config.used.yaml").is_file()),
            }))
    entries.sort(key=lambda e: e[0], reverse=True)
    return jsonify({"runs": [e[1] for e in entries[:50]]})


@app.route("/runs/<job_id>", methods=["DELETE"])
def delete_run(job_id):
    """Delete a finished run's directory. Refuses an in-progress run."""
    job = JOBS.get(job_id)
    if job and not job["done"]:
        return jsonify({"error": "That run is still in progress."}), 409
    d = (CFG["runs_dir"] / job_id).resolve()
    if d.parent != CFG["runs_dir"].resolve() or not d.is_dir():
        return jsonify({"error": "not found"}), 404
    shutil.rmtree(d, ignore_errors=True)
    JOBS.pop(job_id, None)
    return jsonify({"ok": True})


@app.route("/runs", methods=["DELETE"])
def delete_all_runs():
    """Delete every finished run dir. Any run still in progress is skipped."""
    rdir = CFG["runs_dir"]
    with JOB_LOCK:
        live = {j["dir"].name for j in JOBS.values() if not j["done"]}
    deleted = 0
    if rdir.is_dir():
        for d in list(rdir.iterdir()):
            if not d.is_dir() or d.name in live:
                continue
            shutil.rmtree(d, ignore_errors=True)
            JOBS.pop(d.name, None)
            deleted += 1
    return jsonify({"ok": True, "deleted": deleted, "skipped": len(live)})


@app.route("/runs/<job_id>/stop", methods=["POST"])
def stop_run(job_id):
    """Halt the live run by killing its process group. It finalizes as
    'interrupted' (resumable) rather than 'error' -- see _run_job. Idempotent-ish:
    returns 409 if the job isn't currently running."""
    job = JOBS.get(job_id)
    if not job or job["done"]:
        return jsonify({"error": "That run isn't running."}), 409
    job["stopped"] = True
    job["log"].append("[gui] Stop requested -- halting the run.")
    pgid = job.get("pgid")
    if pgid:
        # SIGTERM the whole group (commec + blast/hmmer children); _run_job then
        # observes the exit and marks the run interrupted.
        try:
            os.killpg(pgid, signal.SIGTERM)
        except OSError:
            pass
    return jsonify({"ok": True})


@app.route("/runs/<job_id>/resume", methods=["POST"])
def resume_run(job_id):
    """Continue an interrupted run in place: drop the (possibly truncated)
    output of the step it died on, then re-run commec with -R so the completed
    steps are reused. Refuses if a screen is already running (one at a time)."""
    with JOB_LOCK:
        if any(not j["done"] for j in JOBS.values()):
            return jsonify({"error": "A screen is already running. "
                                     "Wait for it to finish."}), 409
    d = (CFG["runs_dir"] / job_id).resolve()
    if d.parent != CFG["runs_dir"].resolve() or not d.is_dir():
        return jsonify({"error": "not found"}), 404
    meta = _read_meta(d)
    if _effective_status(meta, False) != "interrupted":
        return jsonify({"error": "That run isn't interrupted; nothing to "
                                 "resume."}), 400
    prefix = meta.get("prefix")
    config_path = d / "config.used.yaml"
    if not prefix or not config_path.is_file():
        return jsonify({"error": "Can't resume: the run's config is missing."}), 400
    # The input must live in the run dir so commec writes output here (not next
    # to the input). New runs always have input.fasta in-dir; for a legacy
    # run whose input was external, copy it in.
    indir = d / "input.fasta"
    if indir.is_file():
        fasta = str(indir)
    elif meta.get("fasta") and os.path.isfile(meta["fasta"]):
        try:
            shutil.copy(meta["fasta"], indir)
        except OSError as exc:
            return jsonify({"error": f"Can't resume: {exc}"}), 400
        fasta = str(indir)
    else:
        return jsonify({"error": "Can't resume: the original input is gone."}), 400

    # Redo only the step that was mid-write when killed (commec's -R trusts any
    # non-empty file, so a truncated one would silently corrupt the verdict).
    killed = _drop_interrupted_artifact(d, prefix)
    cmd = _build_command(fasta, d, config_path, prefix, resume=True)
    job = {
        "id": d.name, "label": meta.get("label", d.name), "prefix": prefix,
        "fasta": fasta, "created": meta.get("created", time.time()),
        "finished": None, "cmd": cmd, "dir": d, "log": [],
        "status": "starting", "returncode": None, "done": False,
        "stopped": False, "proc": None, "pgid": None,
    }
    with JOB_LOCK:
        JOBS[d.name] = job
    if killed:
        job["log"].append(f"[gui] Resuming run; re-running the interrupted step "
                          f"(dropped {killed}). Completed steps are reused.")
    else:
        job["log"].append("[gui] Resuming run with -R; existing output reused.")
    job["log"].append("$ " + " ".join(shlex.quote(c) for c in cmd))
    threading.Thread(target=_run_job, args=(job,), daemon=True).start()
    return jsonify({"job_id": d.name})


DESCRIPTION = "Launch the commec screening web GUI (spawns a local Flask server)."


def add_args(parser):
    """Register the GUI server's command-line arguments (subcommand contract)."""
    ap = parser  # local alias keeps the argument definitions below unchanged
    ap.add_argument("--host", default=None,
                    help="Bind address. Default: 127.0.0.1, or 0.0.0.0 if --lan.")
    ap.add_argument("--port", type=int, default=443,
                    help="Port to serve on (default: 443). Binding <1024 needs "
                         "a privilege grant for the unprivileged server -- see "
                         "deploy-notes.md. Use a high port (e.g. 8765) without one.")
    ap.add_argument("--lan", action="store_true",
                    help="Bind 0.0.0.0 so other machines on the LAN can connect.")
    ap.add_argument("--kiosk", action="store_true",
                    help="Auto-open a browser at the server (local kiosk mode).")
    ap.add_argument("--commec-bin", default="commec",
                    help="Path to the commec binary (default: commec on PATH).")
    ap.add_argument("--runs-dir", default="runs",
                    help="Persistent dir holding one directory per run: the "
                         "submitted sequence, commec's search intermediates, "
                         "and the polished results, all kept on disk "
                         "(default: ./runs).")
    ap.add_argument("--runs-keep", type=int, default=0,
                    help="Max retained run dirs (oldest pruned, intermediates "
                         "included). 0 (default) keeps everything.")
    ap.add_argument("--sweep-interval", type=int, default=300,
                    help="Seconds between orphan-process sweeps (default: 300).")
    ap.add_argument("--databases", default="",
                    help="Database directory passed to every screen (via -d).")
    ap.add_argument("--browse-root", default=str(Path.home()),
                    help="Root directory the file browser is confined to "
                         "(default: the user's home directory).")
    ap.add_argument("--password-file",
                    default=str(Path.home() / ".config" / "commec-gui" / "password.hash"),
                    help="File holding the access-password hash (werkzeug "
                         "format). Non-localhost access is blocked until it is "
                         "present + non-empty, then requires login. Default: "
                         "~/.config/commec-gui/password.hash")
    ap.add_argument("--require-local-auth", action="store_true",
                    help="Also require the password from localhost (the walk-up "
                         "kiosk), not just remote clients.")
    ap.add_argument("--threads", type=int, default=0,
                    help="CPU threads passed to commec search tools (-t). "
                         "0 (default) uses all logical cores: one screen runs "
                         "at a time, so the kiosk is dedicated to it.")
    ap.add_argument("--presets", default=str(Path(__file__).parent / "presets.yaml"),
                    help="Path to the presets manifest (default: presets.yaml).")
    ap.add_argument("--tls-cert", default="",
                    help="Path to a TLS certificate (PEM). Enables HTTPS when "
                         "given together with --tls-key.")
    ap.add_argument("--tls-key", default="",
                    help="Path to the TLS private key (PEM) matching --tls-cert.")
    ap.add_argument("--tls-auto", action="store_true",
                    help="Serve HTTPS, generating a per-user cert "
                         "(~/.local/share/commec-gui/certs) on first startup "
                         "(via gen-cert.sh) if it doesn't exist.")
    return parser


def run(args):
    """Configure the server from parsed args and launch it (blocks in app.run)."""
    CFG["commec_bin"] = args.commec_bin
    CFG["runs_dir"] = Path(args.runs_dir).resolve()
    CFG["runs_dir"].mkdir(parents=True, exist_ok=True)
    CFG["runs_keep"] = args.runs_keep
    CFG["default_databases"] = args.databases
    CFG["browse_root"] = Path(args.browse_root).resolve()
    CFG["require_local_auth"] = args.require_local_auth
    try:
        CFG["password_hash"] = (Path(args.password_file)
                                .read_text(encoding="utf-8").strip()) or None
    except OSError:
        CFG["password_hash"] = None
    print("Auth: " + ("ON (password required for "
          + ("all clients" if args.require_local_auth else "non-localhost")
          + ")" if CFG["password_hash"] else
          "LAN blocked (no password file; localhost only)"))
    # 0 => all logical cores (one screen at a time, so dedicate the box to it).
    CFG["threads"] = args.threads or (os.cpu_count() or 1)
    # presets.yaml is deployment-specific (gitignored; the kiosk image provisions
    # it). Fall back to the tracked template so a fresh checkout still runs.
    presets_path = args.presets
    if not os.path.exists(presets_path) and os.path.exists(presets_path + ".example"):
        presets_path += ".example"
        print(f"No {args.presets}; using {presets_path}")
    CFG["presets"] = load_presets(presets_path)
    print(f"Loaded {len(CFG['presets'])} preset(s): "
          + ", ".join(p["id"] for p in CFG["presets"]))
    print(f"commec search tools will use {CFG['threads']} thread(s).")
    print(f"runs dir (kept): {CFG['runs_dir']}"
          + (f", max {CFG['runs_keep']}" if CFG["runs_keep"] else ""))

    # Kill any process group orphaned by a previous (crashed/killed) run, then
    # sweep periodically. Kill in-flight children on shutdown too.
    _sweep()
    threading.Thread(target=_sweep_loop, args=(args.sweep_interval,),
                     daemon=True).start()
    # Auto-mount removable media (read-only) so USB sticks just appear.
    threading.Thread(target=_automount_loop, args=(3,), daemon=True).start()
    atexit.register(_kill_running_jobs)
    for _sig in (signal.SIGTERM, signal.SIGINT):
        signal.signal(_sig, lambda *_a: (_kill_running_jobs(), os._exit(0)))

    # TLS is on only when both cert and key are supplied. Fail loudly on a
    # half-configured or missing pair rather than silently serving cleartext.
    here = Path(__file__).parent
    cert, key = args.tls_cert, args.tls_key

    # --tls-auto: default to a per-user, writable cert dir (XDG data home, i.e.
    # ~/.local/share/commec-gui/certs), generating the pair on first startup if
    # absent. NOT the package dir: a conda/pip install lives in a root-owned or
    # read-only site-packages the server user cannot write. Explicit
    # --tls-cert/--tls-key still take precedence.
    if args.tls_auto and not (cert or key):
        cert_dir = (Path(os.environ.get("XDG_DATA_HOME")
                         or Path.home() / ".local" / "share")
                    / "commec-gui" / "certs")
        cert = str(cert_dir / "server.crt")
        key = str(cert_dir / "server.key")
        if not (os.path.isfile(cert) and os.path.isfile(key)):
            print(f"No TLS cert found; generating one with gen-cert.sh into {cert_dir} ...")
            try:
                subprocess.run(["bash", str(here / "gen-cert.sh")], check=True,
                               env={**os.environ, "COMMEC_GUI_CERT_DIR": str(cert_dir)})
            except (subprocess.CalledProcessError, OSError) as exc:
                raise SystemExit(f"commec gui: error: automatic cert generation "
                                 f"failed: {exc}")

    ssl_context = None
    if cert or key:
        if not (cert and key):
            raise SystemExit("commec gui: error: --tls-cert and --tls-key must "
                             "be given together.")
        for label, path in (("cert", cert), ("key", key)):
            if not os.path.isfile(path):
                raise SystemExit(f"commec gui: error: TLS {label} not found: {path}")
        ssl_context = (cert, key)

    scheme = "https" if ssl_context else "http"
    # Session cookie only over HTTPS when we're serving it. Warn if auth is on
    # but transport is cleartext (login + cookie would be exposed).
    app.config["SESSION_COOKIE_SECURE"] = bool(ssl_context)
    if CFG["password_hash"] and not ssl_context:
        print("WARNING: auth password is set but TLS is OFF -- credentials "
              "would be sent in cleartext. Serve HTTPS in production.")
    host = args.host or ("0.0.0.0" if args.lan else "127.0.0.1")
    # Omit the port from the opened URL when it's the scheme default.
    default_port = 443 if ssl_context else 80
    portpart = "" if args.port == default_port else f":{args.port}"
    url = f"{scheme}://127.0.0.1{portpart}/"
    print(f"commec-gui serving on {scheme}://{host}:{args.port}/  "
          f"(runs dir: {CFG['runs_dir']})")

    if args.kiosk:
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()

    app.run(host=host, port=args.port, threaded=True, ssl_context=ssl_context)


def main():
    """Standalone entrypoint (python -m commec.gui.server / direct script run)."""
    ap = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(ap)
    run(ap.parse_args())


if __name__ == "__main__":
    main()
