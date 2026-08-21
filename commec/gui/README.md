# commec-gui

A tiny web frontend for `commec screen`. One small Flask server serves a single page; the same server powers
both usage modes:

- **Kiosk** (local): binds `127.0.0.1` and opens a browser on this machine.
- **LAN** (remote): binds `0.0.0.0` so other machines on the network connect
  to `https://<this-machine-ip>:<port>/`.

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
commec gui --kiosk --databases /path/to/databases          # local, HTTP
commec gui --lan --tls-auto --databases /path/to/databases # LAN, HTTPS
```

## Configuration (`.env`)

Read by `kiosk.sh` / `lan.sh`. Copy `.env.example` to `.env` (gitignored):

| Variable | Purpose |
|---|---|
| `COMMEC_DB_DIR` | Database directory passed to every screen |
| `COMMEC_CONDA_ENV` | Conda env with commec + flask (default `commec-dev`) |
| `COMMEC_GUI_PORT` | Port to listen on (default `443`; privileged — see TLS/deploy notes. Use a high port like `8765` without a bind grant) |
| `COMMEC_GUI_THREADS` | CPU threads for commec's search tools (`-t`); unset = all cores |
| `COMMEC_GUI_PASSWORD_FILE` | Path to the access-password hash file (default `~/.config/commec-gui/password.hash`); non-localhost access is blocked until it exists, then clients must log in |
| `COMMEC_TLS` | Set to `0` to force plain HTTP in `lan.sh` |
| `COMMEC_TLS_CERT` / `COMMEC_TLS_KEY` | Use a specific cert/key instead of auto-gen |

## Set the LAN password

LAN access stays blocked until a password hash exists. From the `commec/gui`
directory, run:

```bash
./set-password.sh
```

The helper sources `.env`, prompts twice without echoing the password, and writes
the hash to `COMMEC_GUI_PASSWORD_FILE` (or the default path) with mode 0600.
Pass another `.env` file as its first argument. Restart the GUI after changing
the password because the server reads the hash at startup.

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
| `--runs-dir` | `./runs` | Persistent dir holding one directory per run (sequence + intermediates + results), all kept on disk. See [Storage & retention](#storage--retention) |
| `--retention-days` | persisted setting (initially `0`) | Completed-run lifetime in whole days; `0` keeps everything. An explicit value replaces the saved GUI setting |
| `--sweep-interval` | `300` | Seconds between orphan-process sweeps |
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

`presets.yaml` is **gitignored and deployment-specific**.
The repo tracks `presets.yaml.example` as a template, and
the server falls back to it when `presets.yaml` is absent,
so a fresh checkout still runs.
To customise: `cp presets.yaml.example presets.yaml`.

At run time the chosen preset's `config` is written into the run directory as
`config.used.yaml` (a reproducibility record) and passed via
`commec screen -y`. The database directory is injected separately with `-d`
(from `--databases` / `COMMEC_DB_DIR`), so presets stay portable across
machines.

To add or change presets, edit `presets.yaml` (or the `.example`) and restart
the server.

> Caveat: do **not** set `do_cleanup: true` in a preset, as this crashes commec
> at the end of a run (upstream bug in `screen_io.py`). Each job already gets
> its own output directory, so cleanup isn't needed.

## Regulatory jurisdiction

Each run can target one or more regulatory jurisdictions. The default,
**All jurisdictions**, preserves commec's strict default. With specific
countries or region groups selected, applicable control lists use full
compliance and the remaining installed lists use conditional compliance.
Selection does not remove control lists from the screen.
Available jurisdictions reflect the installed control-list data.
Country entries identify their applicable multi-country groups; single-country
aliases are shown once as their ISO country.

The selected codes are stored under
`databases.control_lists.regions` in the run's `config.used.yaml`.

## Inputs

Each run takes a required **run label** (it names the output files and must be
unique) and a sequence supplied one of three ways:

- **Paste**: FASTA text, or rows copied from a spreadsheet (the most DNA-like
  column is taken as the sequence).
- **Server file path**: pick a file via the built-in browser, confined to
  `--browse-root` (plus any read-only USB volume that auto-mounts).
- **Upload**: choose or drag-and-drop a FASTA/CSV/TSV file (200 MB cap).
- **USB key**: Plugging in a flash-drive automatically pops up a tab, from which one can select files (requires some host config).

Every input (paste, upload, spreadsheet, or a picked server file) is normalised
into an `input.fasta` **inside the run dir**. This is required: commec writes its
output next to the input file, so the input must live in the run dir or the
output would scatter into the source directory. Every normalised sequence is
passed to `commec screen`, regardless of length.

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

### Interrupted runs (resume)

If the **server process** dies mid-run (crash, reboot, kill), that run's marker
is stuck at `running`; on the next startup the orphan sweep promotes it to
**Interrupted** (and reaps any still-running search process). Interrupted runs
show an **Interrupted** badge and a **Resume** button. Resuming re-runs commec
with `-R` in the same directory, re-using the completed search steps, but
first it **deletes the newest step output** (`output_<prefix>/` or
`input_<prefix>/`), because that's the step that could have been mid-write when the process
died and commec's `-R` would otherwise trust a truncated file (a silent way to
lose a hit). So the cut-off step is recomputed while everything before it is
reused. (A run that commec itself *errors* on, or that the user *cancels*,
finalizes as `error`, not `interrupted`: those aren't offered for resume.)

## Storage & retention

Each run gets one persistent directory under `--runs-dir` (named by job id) that
holds **everything**: the submitted sequence (`input.fasta`), `config.used.yaml`,
commec's raw search intermediates (`output_<prefix>/`, `input_<prefix>/`), and the
polished artifacts (`<prefix>.output.json`, `<prefix>_summary.html`,
`<prefix>.screen.log`).

Raw query sequences persist on disk. That is a deliberate
trade-off: the kiosk box is treated as secured, so on-box debuggability is worth
more than ephemerality (see `deploy-notes.md`). The GUI's Results panel shows the
polished artifacts by default and reveals the intermediates (and `input.fasta`)
under the **Advanced options** toggle.

Retention is forever by default. The Results card stores a server-wide retention
period in `--runs-dir/.retention.json`; `0` keeps everything. A terminal run
expires when its recorded completion time is more than that many whole days old.
Expiry removes the whole run directory, including inputs and intermediates. Live
and nonterminal runs are excluded. `--retention-days N` can set and persist the
same value at startup.

A periodic sweep (`--sweep-interval`) reaps process groups orphaned by a
crash/restart, flags their runs **Interrupted** (resumable, see above), and
expires completed runs.

## Planned

Soft launch-into-fullscreen for the kiosk (`COMMEC_GUI_FULLSCREEN`, see
`deploy-notes.md`); the in-app Fullscreen button already works.
