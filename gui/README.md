# commec-gui

A tiny web frontend for `commec screen` One small Flask server serves a single page; the same server powers
both usage modes:

- **Kiosk** (local): binds `127.0.0.1` and opens a browser on this machine.
- **LAN** (remote): binds `0.0.0.0` so other machines on the network connect
  to `http://<this-machine-ip>:<port>/`.

The server shells out to the `commec` binary, so run it inside the
`commec-dev` conda env (or point `--commec-bin` at the binary).

## Run

### Kiosk mode

Copy `.env.example` to `.env`, set `COMMEC_DB_DIR` (and optionally the conda
env name and port), then launch:

```bash
cp .env.example .env   # then edit COMMEC_DB_DIR
./kiosk.sh
```

`kiosk.sh` sources `.env`, activates the conda env, and starts the server in
kiosk mode (localhost + auto-open browser) with the database directory
prefilled in the form. Pass a different env file as the first argument:
`./kiosk.sh /path/to/other.env`.

### LAN mode (HTTPS)

```bash
./lan.sh
```

`lan.sh` sources `.env`, binds `0.0.0.0`, and serves **HTTPS** by default,
auto-generating `certs/server.{crt,key}` on first startup (see TLS below).
Clients connect at `https://<this-machine-ip>:<port>/`.

### Manual

```bash
conda activate commec-dev
python server.py --kiosk --databases /path/to/databases          # local, HTTP
python server.py --lan --tls-auto --databases /path/to/databases # LAN, HTTPS
```

## Configuration (`.env`)

Read by `kiosk.sh` / `lan.sh`. Copy `.env.example` to `.env` (gitignored):

| Variable | Purpose |
|---|---|
| `COMMEC_DB_DIR` | Database directory passed to every screen |
| `COMMEC_CONDA_ENV` | Conda env with commec + flask (default `commec-dev`) |
| `COMMEC_GUI_PORT` | Port to listen on (default `443`; privileged — see TLS/deploy notes. Use a high port like `8765` without a bind grant) |
| `COMMEC_GUI_THREADS` | CPU threads for commec's search tools (`-t`); unset = all cores |
| `COMMEC_GUI_PASSWORD_FILE` | Path to the access-password hash file (default `~/.config/commec-gui/password.hash`); if present, non-localhost clients must log in |
| `COMMEC_TLS` | Set to `0` to force plain HTTP in `lan.sh` |
| `COMMEC_TLS_CERT` / `COMMEC_TLS_KEY` | Use a specific cert/key instead of auto-gen |

## TLS / encryption

LAN traffic (including the FASTA sequences) is encrypted over HTTPS. TLS only
encrypts the transport; it does **not** authenticate clients. Anyone on the
LAN can still reach the GUI.

- `lan.sh` (or `--tls-auto`) generates `certs/server.{crt,key}` on first
  startup via `gen-cert.sh`, which **requires mkcert**. The cert's SANs cover
  localhost, the hostname, and the detected LAN IP. If mkcert isn't installed,
  `gen-cert.sh` errors out (set `COMMEC_TLS=0` to serve plain HTTP for
  local/dev instead).
- **To remove the browser "not trusted" warning in kiosk mode**: run `mkcert -install`. It adds the local CA to the trust stores. It needs`libnss3-tools` for Firefox's store, and sudo for the system store. For other
  machines on the LAN, you can get rid of the warning by installing the root CA printed by `gen-cert.sh`
  (`"$(mkcert -CAROOT)/rootCA.pem"`) into each client's trust store.


## Server options

| Flag | Default | Purpose |
|---|---|---|
| `--host` | `127.0.0.1` (or `0.0.0.0` with `--lan`) | Bind address |
| `--port` | `443` | Port (privileged; needs a bind grant for the unprivileged server, or use a high port) |
| `--lan` | off | Bind `0.0.0.0` for LAN access |
| `--kiosk` | off | Auto-open a local browser |
| `--tls-auto` | off | Serve HTTPS, generating a cert on first startup if absent |
| `--tls-cert` / `--tls-key` | (none) | Use a specific cert/key pair (HTTPS) |
| `--commec-bin` | `commec` | Path to the commec binary |
| `--work-dir` | `./runs` | Scratch for live runs (submitted sequences + intermediates); deleted at end-of-run. Put it on tmpfs in prod. See [Storage & retention](#storage--retention) |
| `--results-dir` | `./results` | Persistent dir for retained results (no raw sequence). Must differ from `--work-dir` |
| `--results-keep` | `0` (unlimited) | Max retained result dirs; oldest pruned past the cap |
| `--sweep-interval` | `300` | Seconds between orphan-scratch sweeps |
| `--databases` | (none) | Database dir passed to every screen (via `-d`) |
| `--browse-root` | home dir | Root the server-side file browser is confined to |
| `--threads` | all logical cores | CPU threads passed to commec search tools (via `-t`) |
| `--presets` | `presets.yaml` | Path to the presets manifest |

## Presets

The GUI is preset-driven: the user supplies a sequence and picks one screening
**preset**. Presets are defined in `presets.yaml`, each with a required `id`
and `label`, an optional `recommended` flag (the GUI shows a recommended /
NOT-RECOMMENDED badge), an optional `description`, and a `config` block using
commec's own config keys (`skip_taxonomy_search`, `skip_nt_search`,
`protein_search_tool`, etc.).

At run time the chosen preset's `config` is written into the run directory as
`config.used.yaml` (a reproducibility record) and passed via
`commec screen -y`. The database directory is injected separately with `-d`
(from `--databases` / `COMMEC_DB_DIR`), so presets stay portable across
machines.

To add or change presets, edit `presets.yaml` and restart the server.

> Caveat: do **not** set `do_cleanup: true` in a preset, as this crashes commec
> at the end of a run (upstream bug in `screen_io.py`). Each job already gets
> its own output directory, so cleanup isn't needed.

## Inputs

Each run takes a required **run label** (it names the output files and must be
unique) and a sequence supplied one of three ways:

- **Paste**: FASTA text, or rows copied from a spreadsheet (the most DNA-like
  column is taken as the sequence).
- **Server file path**: pick a file via the built-in browser, confined to
  `--browse-root` (plus any read-only USB volume that auto-mounts).
- **Upload**: choose or drag-and-drop a FASTA/CSV/TSV file (200 MB cap).
- **USB key**: Plugging in a flash-drive automatically pops up a tab, from which one can select files (requires some host config).

Pasted/uploaded/spreadsheet input is normalised into a FASTA file in the run
dir; a server-side file path is screened in place. A default-on **"skip
sequences shorter than 50 bp"** toggle drops sub-threshold records.

One screen runs at a time (screening is CPU-bound and `--threads` already
saturates the box). The shared run state is polled by every connected browser,
so a second client reflects the busy/idle status and attaches to an in-progress
run's live log.

## Results

When a run finishes the page shows a **verdict summary** (per-query
Flag/Warning/Pass/Skip, parsed from commec's `*.output.json`), a link to open
the embedded HTML **report**, and the list of retained **output files**. The
**Results** panel lists finished runs (newest first) and survives a server
restart; each can be reopened or deleted.

## Storage & retention

Sequences are sensitive; results are not, so the two live in different places:

- **Work dir (`--work-dir`, sensitive + ephemeral).** Holds the submitted
  sequences and commec's raw search intermediates. Wiped at end-of-run, with a
  periodic sweep (`--sweep-interval`) reaping orphans from crashed/killed runs.
  **Put it on tmpfs in production** so raw sequence never touches persistent
  disk (see `deploy-notes.md`).
- **Results dir (`--results-dir`, retained).** Only the polished, self-
  contained artifacts are copied off the work dir and kept:
  `<prefix>.output.json`, `<prefix>_summary.html`, `<prefix>.screen.log`, and
  `config.used.yaml` -- none contain raw query sequence. Retention is forever by
  default; set `--results-keep N` for a rolling cap.

The two directories must differ (the server refuses to start otherwise).

## Planned

Soft launch-into-fullscreen for the kiosk (`COMMEC_GUI_FULLSCREEN`, see
`deploy-notes.md`); the in-app Fullscreen button already works.
