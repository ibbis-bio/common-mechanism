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

import yaml
from flask import (Flask, Response, jsonify, redirect, request, send_file,
                   session, url_for)
from werkzeug.security import check_password_hash

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
    # work_dir holds the live run scratch (sensitive sequences + intermediates);
    # put it on tmpfs in prod. results_dir holds the retained polished results
    # (no raw sequence) and should be persistent disk.
    "work_dir": Path("runs").resolve(),
    "results_dir": Path("results").resolve(),
    "results_keep": 0,      # max retained result dirs; 0 = unlimited (forever)
    "default_databases": "",
    "threads": None,        # CPU threads passed to commec (-t); None = don't set
    "browse_root": Path.home(),  # file browser is confined to this directory
    "password_hash": None,  # set from --password-file; None = no auth
    "require_local_auth": False,  # also require the password from localhost
    "presets": [],          # list of {id, label, description, config}
}

# Polished result artifacts we persist off the (tmpfs) scratch. Everything else
# in a run dir is sensitive (query sequences) or bulky raw search output, and is
# deleted at end-of-run. {prefix} is filled with the run's output prefix.
RESULT_ARTIFACTS = (
    "{prefix}.output.json",
    "{prefix}_summary.html",
    "{prefix}.screen.log",
    "config.used.yaml",
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


# Top-level database subdirectories under the database dir, with friendly
# labels for messages. Mirrors commec's screen-default-config.yaml layout.
DB_LABELS = {
    "biorisk": "biorisk",
    "low_concern": "low-concern",
    "nr_blast": "regulated-protein (blastx)",
    "nr_dmnd": "regulated-protein (diamond)",
    "nt_blast": "regulated-nucleotide",
    "taxonomy": "taxonomy",
}


def _required_db_dirs(config):
    """Which db subdirectories a preset needs, derived from its skip flags.

    Mirrors commec's gating: protein search runs unless taxonomy is skipped;
    nucleotide search runs only when neither taxonomy nor nt is skipped.
    """
    dirs = ["biorisk", "low_concern"]
    skip_tx = bool(config.get("skip_taxonomy_search", False))
    skip_nt = bool(config.get("skip_nt_search", False))
    tool = config.get("protein_search_tool", "blastx")
    if not skip_tx:
        dirs.append("nr_dmnd" if tool == "diamond" else "nr_blast")
        dirs.append("taxonomy")
        if not skip_nt:
            dirs.append("nt_blast")
    return dirs


def _missing_databases(preset):
    """Return the required db subdirs that are absent under the database dir.

    Only meaningful when a database dir is configured and the preset relies on
    it (our presets don't declare their own db paths). Returns [] otherwise.
    """
    db_dir = CFG["default_databases"]
    if not db_dir or not os.path.isdir(db_dir):
        return []
    return [name for name in _required_db_dirs(preset["config"])
            if not os.path.isdir(os.path.join(db_dir, name))]

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
        return None
    total = info.get("MemTotal")
    avail = info.get("MemAvailable", info.get("MemFree", 0))
    if not total:
        return None
    used = total - avail
    return (100.0 * used / total, used, total)


@app.route("/metrics")
def metrics():
    cpu = _cpu_percent()
    mem = _mem_info()
    return jsonify({
        "cpu": cpu,
        "ram": ({"pct": mem[0], "used": mem[1], "total": mem[2]} if mem else None),
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
# into a FASTA file in the run dir, optionally dropping sub-threshold sequences.
# ---------------------------------------------------------------------------
MIN_SEQ_LEN = 50  # "skip short sequences" threshold (bp), default-on in the UI


def _clean_seq(s):
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
            seq.append(_clean_seq(line))
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
        records.append((names[0] if names else f"sequence_{i + 1}", _clean_seq(seq)))
    return records


def _parse_blocks_text(text):
    """Header-less raw sequences: blocks separated by a blank line, each block's
    lines joined into one sequence (so a single sequence may wrap over lines)."""
    records = []
    for i, block in enumerate(re.split(r"\n\s*\n", text.strip())):
        seq = _clean_seq(block)
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
        seq = _clean_seq(seq)
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


def _build_command(preset, fasta, outdir, config_path, prefix):
    """Build the commec screen argv list for a preset run.

    Args are always passed as a list (never shell=True) so user-supplied
    paths can never be interpreted as shell syntax. The preset's config is
    passed via -y; the database directory is injected via -d so it overrides
    whatever (if anything) the preset declares.
    """
    cmd = [CFG["commec_bin"], "screen", "-y", str(config_path)]
    if CFG["default_databases"]:
        cmd += ["-d", CFG["default_databases"]]
    if CFG["threads"]:
        cmd += ["-t", str(CFG["threads"])]
    # -o as a basename (outdir/prefix, no trailing sep) names the outputs after
    # the run label, independent of the input filename. Each job gets its own
    # fresh output directory, so no force/resume needed.
    cmd += ["-o", str(outdir / prefix), fasta]
    return cmd


def _finalize_job(job):
    """Persist the polished results off the scratch dir, then delete the
    scratch dir entirely (removing the sensitive sequences + raw intermediates).

    Writes meta.json (status/label/verdict) as a completion marker that the
    sweep uses to tell a clean finish from an orphan, and that survives a Flask
    restart for the results view.
    """
    scratch, resultdir, prefix = job["scratch"], job["resultdir"], job["prefix"]
    try:
        resultdir.mkdir(parents=True, exist_ok=True)
        for tmpl in RESULT_ARTIFACTS:
            src = scratch / tmpl.format(prefix=prefix)
            if src.is_file():
                shutil.copy2(src, resultdir / src.name)
        meta = {
            "id": job["id"], "label": job["label"], "prefix": prefix,
            "status": job["status"], "returncode": job["returncode"],
            "created": job["created"], "finished": time.time(),
            "summary": _summarize_output(resultdir),
        }
        (resultdir / "meta.json").write_text(json.dumps(meta), encoding="utf-8")
    except OSError as exc:
        job["log"].append(f"[gui] WARNING: could not persist results: {exc}")
    # Delete the sensitive scratch (inputs, translations, raw search output).
    shutil.rmtree(scratch, ignore_errors=True)
    job["finished"] = time.time()
    _sweep_results()  # apply any rolling-retention cap


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
        (job["scratch"] / ".pgid").write_text(str(proc.pid), encoding="utf-8")
    except OSError:
        pass
    job["status"] = "running"
    for line in proc.stdout:
        job["log"].append(line.rstrip("\n"))
    proc.wait()
    job["returncode"] = proc.returncode
    job["status"] = "done" if proc.returncode == 0 else "error"
    # Persist results + wipe scratch BEFORE marking done, so a client that loads
    # results on the 'end' event reads the persisted copy.
    _finalize_job(job)
    job["done"] = True


# ---------------------------------------------------------------------------
# Cleanup: orphan sweep (crashed/killed runs) and process-group killing.
# ---------------------------------------------------------------------------
def _kill_pgid(pgid):
    try:
        os.killpg(int(pgid), signal.SIGKILL)
    except (ProcessLookupError, ValueError, OSError):
        pass


def _sweep_scratch():
    """Delete scratch dirs not owned by a live run (orphans from crashes or a
    server restart), killing any still-running process group first."""
    work = CFG["work_dir"]
    if not work.is_dir():
        return
    with JOB_LOCK:
        live = {j["scratch"].name for j in JOBS.values() if not j["done"]}
    for d in work.iterdir():
        if not d.is_dir() or d.name in live:
            continue
        pgid_file = d / ".pgid"
        if pgid_file.is_file():
            try:
                _kill_pgid(pgid_file.read_text(encoding="utf-8").strip())
            except OSError:
                pass
        shutil.rmtree(d, ignore_errors=True)


def _sweep_results():
    """Apply the rolling-retention cap, if any (default: keep everything)."""
    keep = CFG["results_keep"]
    if not keep or keep <= 0 or not CFG["results_dir"].is_dir():
        return
    dirs = [d for d in CFG["results_dir"].iterdir() if d.is_dir()]
    dirs.sort(key=lambda d: d.stat().st_mtime, reverse=True)
    for d in dirs[keep:]:
        shutil.rmtree(d, ignore_errors=True)


def _sweep():
    _sweep_scratch()
    _sweep_results()


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
# Authentication. When a password hash is configured, every NON-localhost
# request needs a valid session. The walk-up kiosk (127.0.0.1) is exempt --
# physical presence is the trust -- unless --require-local-auth is set.
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
    """Expose the preset choices the page renders, with DB availability."""
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
    return jsonify({"presets": presets})


def _label_exists(label):
    """True if a retained run already uses this label (keeps labels unique)."""
    rdir = CFG["results_dir"]
    if not rdir.is_dir():
        return False
    for d in rdir.iterdir():
        if d.is_dir() and _read_meta(d).get("label") == label:
            return True
    return False


@app.route("/screen", methods=["POST"])
def screen():
    preset_id = (request.form.get("preset") or "").strip()
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

    # Resolve the input into either a pass-through path (server file, unchanged)
    # or a set of records we materialise into the run dir.
    mode = (request.form.get("input_mode") or "path").strip()
    skip_short = (request.form.get("skip_short") or "1") != "0"
    records, passthrough = None, None

    if mode == "paste":
        records = parse_pasted_sequences(request.form.get("sequence_text") or "")
        if not records:
            return jsonify({"error": "Couldn't find any sequences in the pasted "
                "text. Paste FASTA, or spreadsheet rows with a sequence column."}), 400
    elif mode == "upload":
        f = request.files.get("upload_file")
        if not f or not f.filename:
            return jsonify({"error": "Choose a file to upload."}), 400
        text = f.read().decode("utf-8", errors="replace")
        records = parse_pasted_sequences(text)
        if not records:
            return jsonify({"error": "No sequences found in the uploaded file. "
                "Upload a FASTA file (or a CSV/TSV with a sequence column)."}), 400
    elif mode == "path":
        fasta = (request.form.get("fasta_file") or "").strip()
        if not fasta:
            return jsonify({"error": "A FASTA file path is required."}), 400
        # Confine to the same roots the browser exposes (home + USB mounts), so
        # a hand-crafted POST can't screen arbitrary server files. Resolve first
        # to defeat symlink/.. escapes, matching /browse and /results/file.
        target = Path(fasta).resolve()
        if not any(target == r or r in target.parents for r in _allowed_roots()):
            return jsonify({"error": "That file is outside the allowed directory."}), 400
        fasta = str(target)
        if not os.path.isfile(fasta):
            return jsonify({"error": f"FASTA path not found on server: {fasta}"}), 400
        if skip_short:
            try:
                records = _parse_fasta_text(Path(fasta).read_text(encoding="utf-8"))
            except OSError as exc:
                return jsonify({"error": f"Could not read {fasta}: {exc}"}), 400
            if not records:
                return jsonify({"error": f"No FASTA records found in {fasta}."}), 400
        else:
            passthrough = fasta
    else:
        return jsonify({"error": f"Unknown input mode: {mode!r}"}), 400

    note = []
    if records is not None:
        records = _normalise_records(records)
        if skip_short:
            kept = [r for r in records if len(r[1]) >= MIN_SEQ_LEN]
            dropped = len(records) - len(kept)
            records = kept
            if dropped:
                note.append(f"Skipped {dropped} sequence(s) shorter than {MIN_SEQ_LEN} bp.")
        if not records:
            return jsonify({"error":
                f"No sequences left to screen (all shorter than {MIN_SEQ_LEN} bp)."}), 400

    with JOB_LOCK:
        busy = [j for j in JOBS.values() if not j["done"]]
        if busy:
            return jsonify({"error": "A screen is already running. "
                                     "Wait for it to finish."}), 409

        job_id = uuid.uuid4().hex[:12]
        scratch = CFG["work_dir"] / job_id      # tmpfs: sensitive, ephemeral
        resultdir = CFG["results_dir"] / job_id  # persistent: retained results
        scratch.mkdir(parents=True, exist_ok=True)

        # Write the preset's config into the scratch dir (also a reproducibility
        # record, persisted with the results) and point commec at it.
        config_path = scratch / "config.used.yaml"
        with open(config_path, "w", encoding="utf-8") as fh:
            yaml.safe_dump(preset["config"], fh, default_flow_style=False)

        if records is not None:
            fasta = str(scratch / "input.fasta")
            _write_fasta(records, fasta)
            note.insert(0, f"Prepared {len(records)} sequence(s) for screening.")
        else:
            fasta = passthrough

        cmd = _build_command(preset, fasta, scratch, config_path, prefix)
        job = {
            "id": job_id,
            "label": label or prefix,
            "prefix": prefix,
            "created": time.time(),
            "finished": None,
            "cmd": cmd,
            "scratch": scratch,
            "resultdir": resultdir,
            "log": [],
            "status": "starting",
            "returncode": None,
            "done": False,
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
    sev = {"Flag": 3, "Warning": 2, "Pass": 1, "Skip": 0}
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
    overall = max((r["status"] for r in rows), key=lambda s: sev.get(s, 0))
    return {"n": len(rows), "overall": overall, "queries": rows}


def _job_dir(job_id):
    """The directory currently holding a run's files: the live (tmpfs) scratch
    while running, the persisted results dir once finished, or -- when the job
    isn't in memory (e.g. after a restart) -- the results dir on disk if it
    exists. Returns None if nothing is found."""
    job = JOBS.get(job_id)
    if job:
        d = job["resultdir"] if job["done"] else job["scratch"]
        return d if d.is_dir() else None
    d = CFG["results_dir"] / job_id
    return d if d.is_dir() else None


@app.route("/results/<job_id>")
def results(job_id):
    job = JOBS.get(job_id)
    d = _job_dir(job_id)
    if d is None and not job:
        return jsonify({"error": "unknown job"}), 404
    files = []
    if d is not None:
        for p in sorted(d.iterdir()):
            if p.is_file() and p.name != "meta.json":
                files.append({"name": p.name, "size": p.stat().st_size})
    if job:
        status, returncode, label = job["status"], job["returncode"], job["label"]
    else:
        meta = _read_meta(d)
        status = meta.get("status", "done")
        returncode = meta.get("returncode")
        label = meta.get("label")
    return jsonify({
        "status": status,
        "returncode": returncode,
        "label": label,
        "files": files,
        "summary": _summarize_output(d) if d is not None else None,
    })


def _read_meta(d):
    try:
        return json.loads((d / "meta.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


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
    persisted meta.json markers so the list survives restarts."""
    rdir = CFG["results_dir"]
    entries = []
    if rdir.is_dir():
        for d in rdir.iterdir():
            if not d.is_dir():
                continue
            meta = _read_meta(d)
            finished = meta.get("finished") or d.stat().st_mtime
            summary = meta.get("summary") or {}
            entries.append((finished, {
                "id": meta.get("id", d.name),
                "label": meta.get("label", d.name),
                "status": meta.get("status", "done"),
                "returncode": meta.get("returncode"),
                "finished": finished,
                "overall": summary.get("overall"),
                "n": summary.get("n"),
            }))
    entries.sort(key=lambda e: e[0], reverse=True)
    return jsonify({"runs": [e[1] for e in entries[:50]]})


@app.route("/runs/<job_id>", methods=["DELETE"])
def delete_run(job_id):
    """Delete a finished run's retained results. Refuses an in-progress run."""
    job = JOBS.get(job_id)
    if job and not job["done"]:
        return jsonify({"error": "That run is still in progress."}), 409
    d = (CFG["results_dir"] / job_id).resolve()
    if d.parent != CFG["results_dir"].resolve() or not d.is_dir():
        return jsonify({"error": "not found"}), 404
    shutil.rmtree(d, ignore_errors=True)
    JOBS.pop(job_id, None)
    return jsonify({"ok": True})


def main():
    ap = argparse.ArgumentParser(description="Tiny web GUI for `commec screen`.")
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
    ap.add_argument("--work-dir", default="runs",
                    help="Scratch dir for live runs; holds the submitted "
                         "sequences and intermediates, deleted at end-of-run. "
                         "Put it on tmpfs in prod (default: ./runs).")
    ap.add_argument("--results-dir", default="results",
                    help="Persistent dir for retained results (no raw sequence) "
                         "(default: ./results).")
    ap.add_argument("--results-keep", type=int, default=0,
                    help="Max retained result dirs (oldest pruned). "
                         "0 (default) keeps everything.")
    ap.add_argument("--sweep-interval", type=int, default=300,
                    help="Seconds between orphan sweeps (default: 300).")
    ap.add_argument("--databases", default="",
                    help="Database directory passed to every screen (via -d).")
    ap.add_argument("--browse-root", default=str(Path.home()),
                    help="Root directory the file browser is confined to "
                         "(default: the user's home directory).")
    ap.add_argument("--password-file",
                    default=str(Path.home() / ".config" / "commec-gui" / "password.hash"),
                    help="File holding the access-password hash (werkzeug "
                         "format). If present + non-empty, non-localhost "
                         "requests must log in. Default: "
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
                    help="Serve HTTPS, generating certs/server.{crt,key} on "
                         "first startup (via gen-cert.sh) if they don't exist.")
    args = ap.parse_args()

    CFG["commec_bin"] = args.commec_bin
    CFG["work_dir"] = Path(args.work_dir).resolve()
    CFG["results_dir"] = Path(args.results_dir).resolve()
    if CFG["work_dir"] == CFG["results_dir"]:
        ap.error("--work-dir and --results-dir must differ (sensitive scratch "
                 "vs retained results).")
    CFG["work_dir"].mkdir(parents=True, exist_ok=True)
    CFG["results_dir"].mkdir(parents=True, exist_ok=True)
    CFG["results_keep"] = args.results_keep
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
          + ")" if CFG["password_hash"] else "off (no password file)"))
    # 0 => all logical cores (one screen at a time, so dedicate the box to it).
    CFG["threads"] = args.threads or (os.cpu_count() or 1)
    CFG["presets"] = load_presets(args.presets)
    print(f"Loaded {len(CFG['presets'])} preset(s): "
          + ", ".join(p["id"] for p in CFG["presets"]))
    print(f"commec search tools will use {CFG['threads']} thread(s).")
    print(f"work dir (scratch): {CFG['work_dir']}")
    print(f"results dir (kept): {CFG['results_dir']}"
          + (f", max {CFG['results_keep']}" if CFG["results_keep"] else ""))

    # Clear any orphan scratch from a previous (crashed/killed) run, then sweep
    # periodically. Kill in-flight children on shutdown so we don't orphan them.
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

    # --tls-auto: default to certs/server.{crt,key}, generating them on first
    # startup if absent. Explicit --tls-cert/--tls-key still take precedence.
    if args.tls_auto and not (cert or key):
        cert = str(here / "certs" / "server.crt")
        key = str(here / "certs" / "server.key")
        if not (os.path.isfile(cert) and os.path.isfile(key)):
            print("No TLS cert found; generating one with gen-cert.sh ...")
            try:
                subprocess.run(["bash", str(here / "gen-cert.sh")], check=True)
            except (subprocess.CalledProcessError, OSError) as exc:
                ap.error(f"Automatic cert generation failed: {exc}")

    ssl_context = None
    if cert or key:
        if not (cert and key):
            ap.error("--tls-cert and --tls-key must be given together.")
        for label, path in (("cert", cert), ("key", key)):
            if not os.path.isfile(path):
                ap.error(f"TLS {label} not found: {path}")
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
          f"(work dir: {CFG['work_dir']})")

    if args.kiosk:
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()

    app.run(host=host, port=args.port, threaded=True, ssl_context=ssl_context)


if __name__ == "__main__":
    main()
