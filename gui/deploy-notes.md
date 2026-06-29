# Deployment / image-minting notes

Notes for the image-minting pipeline that will provision the kiosk image and
write the GUI's `.env`. Keep this in sync with `server.py` / `kiosk.sh` /
`lan.sh`.

## Runtime environment requirement (important)

The GUI server (`server.py`) does **not** import commec; it shells out to the
`commec` binary. But commec only works when its whole toolchain resolves on
`PATH`, so the server process must run in an environment that has **all** of:

- `commec`
- search tools: `blastx`, `blastn`, `diamond`, `hmmscan`, `cmscan`, `transeq`,
  `makeblastdb`
- `taxonkit` -- commec imports `pytaxonkit`, which runs `taxonkit version` at
  **import time**. If `taxonkit` is not on `PATH`, commec dies before doing
  anything (FileNotFoundError), even for the biorisk-only preset.

`kiosk.sh` / `lan.sh` achieve this by `conda activate "$COMMEC_CONDA_ENV"`
before launching, which puts that env's `bin/` on `PATH`. The image must point
`COMMEC_CONDA_ENV` at the env that actually contains commec + the tools above.

The server itself additionally needs `flask` and `pyyaml` installed in that
same env.

## `.env` values the image should set

`kiosk.sh` / `lan.sh` source `.env`. The image should write:

| Variable | Value on the Debian prototype | Notes |
|---|---|---|
| `COMMEC_CONDA_ENV` | `commec-env` | NOT the script default `commec-dev`. Must be the env with commec + tools + taxonkit + flask + pyyaml. |
| `COMMEC_DB_DIR` | `/home/<user>/commec-dbs` | Database dir. Must contain `biorisk/`, `low_concern/`, `nr_blast/` (or `nr_dmnd/`), `nt_blast/`, `taxonomy/` with the filenames commec's `screen-default-config.yaml` expects. |
| `COMMEC_GUI_PORT` | `8765` | |
| `COMMEC_GUI_THREADS` | unset | CPU threads for commec's search tools (`-t`). Unset = all logical cores. Only set to cap usage (e.g. shared host). |
| `COMMEC_GUI_FULLSCREEN` | `1` (kiosk) | Planned (see below). When set, the page auto-enters fullscreen on first interaction. |
| TLS | (see below) | |

## Fullscreen display

Two pieces, by design:

- **In-app toggle button** (DOM Fullscreen API) -- a Fullscreen/Exit button in
  the page. Works everywhere; reversible. No image config needed.
  **Implemented** (header button in `index.html`).
- **"Soft" launch-into-fullscreen for the kiosk** (planned -- not yet wired) --
  this is the part the image pipeline configures:
  - Set `COMMEC_GUI_FULLSCREEN=1` in `.env`. The server surfaces this to the
    page, which then calls `requestFullscreen()` on the **first** user gesture
    (the operator's first touch/click). The browser blocks auto-fullscreen on
    load without a gesture, so this is the robust way to "boot into fullscreen"
    while keeping the toggle button functional.
  - Launch the browser **maximized/normal -- do NOT use `firefox-esr --kiosk`**.
    `--kiosk` is a window-level lock that the in-app Exit button can't undo;
    it's only for a fully locked-down public terminal (a separate, "hard"
    option we're not using here).

## Privileged port (443)

The GUI now defaults to **port 443** so clients get a clean `https://<host>/`
URL with no port to type. 443 is a privileged port (<1024), and the server
runs as the unprivileged GUI user, so the host must grant the bind. Pick one
(the image controls this):

- **systemd service unit (recommended):** run the GUI as a service with
  `User=<kiosk-user>` and `AmbientCapabilities=CAP_NET_BIND_SERVICE`. Grants
  exactly the bind to just this service; nothing else gains privilege.
- **Lower the unprivileged-port floor (simple, dedicated box):**
  `sysctl -w net.ipv4.ip_unprivileged_port_start=443` (persist in
  `/etc/sysctl.d/`). Any process may then bind 443+. Fine for a single-purpose
  appliance.
- **setcap on the interpreter (blunter):**
  `setcap cap_net_bind_service=+ep "$(readlink -f .../bin/python3)"` -- grants
  the cap to that interpreter for all scripts it runs.

Without any grant, set `COMMEC_GUI_PORT`/`--port` to a high port (e.g. 8765).
Optional nicety: redirect 80 -> 443 (nft/iptables) so the bare hostname
upgrades to HTTPS.

## Removable media (USB) auto-mount

The GUI auto-detects removable drives so users get the Windows/Mac "plug it in
and it's just there" behaviour. A background thread polls `lsblk` for removable
partitions (`RM=1`) and mounts any unmounted one **read-only** via
`udisksctl mount -o ro`; mounted volumes then appear as a tab in the page. The
server only ever READS the media -- it never writes to or executes from it.

`udisksctl mount` must be **authorized for the server process**. On the
prototype (server launched over ssh, outside the graphical session) polkit
denies it by default (`NotAuthorizedCanObtain`). The packaging workflow must
grant it, by EITHER:

- **Launching the GUI from the kiosk's graphical session** (autostart in the
  logged-in seat) -- udisks2 mounts are then authorized for the active local
  session with no extra rule (the GUI user should be in `plugdev`); OR
- **Installing a polkit rule** allowing `plugdev` to mount removable media,
  which works regardless of how the server is launched:

  `/etc/polkit-1/rules.d/50-commec-usb.rules`
  ```
  polkit.addRule(function(action, subject) {
    if (action.id.indexOf("org.freedesktop.udisks2.filesystem-mount") == 0 &&
        subject.isInGroup("plugdev")) { return polkit.Result.YES; }
  });
  ```

Mounts are read-only. For extra hardening against hostile filesystems, also
apply `nodev,nosuid,noexec` to removable mounts (udisks2/mount config). If you
prefer, a system-level auto-mount (udev+systemd / usbmount) works too -- then
the server just reads `/media` and the server-side mount step is a no-op.

## TLS for production

LAN mode serves HTTPS. `--tls-auto` (via `gen-cert.sh`) issues a **per-kiosk**
cert generated **on the device** at first boot -- so **no private key is ever
baked into the image** (we assume the image leaks). It needs no network, so it
works on fully-offline kiosks.

Decision: we are **not** using certbot / a public CA. Offline kiosks can't do
ACME issuance or 90-day renewals, and DNS-01 credentials would be a baked
secret -- both disqualifying. The chosen approach:

- Install `mkcert` (+ `libnss3-tools` so Firefox's NSS store is updated) in the
  image, and at first boot run `mkcert -install` as the kiosk user. `gen-cert.sh`
  then issues an mkcert-signed cert. The kiosk's **own** browser trusts it ->
  no warning on the walk-up screen, fully automated, offline, nothing baked.
- **Remote LAN clients (laptops) still get a warning** until that kiosk's root
  CA (`$(mkcert -CAROOT)/rootCA.pem`) is installed in their trust store. That's
  a local-IT task we deliberately do not automate from the kiosk; if a site
  dislikes the warning, their IT installs the CA.
- Without mkcert, `gen-cert.sh` **errors out** (set `COMMEC_TLS=0` for plain
  HTTP) -- there is no self-signed openssl fallback, so a missing mkcert is a
  loud provisioning failure rather than a silently-shipped warning-cert. The
  image must also install `libnss3-tools` so `mkcert -install` can populate
  Firefox's NSS trust store (otherwise the cert is mkcert-signed but the
  kiosk's own Firefox still warns).

Each kiosk having its own CA means a compromised device exposes only its own
CA -- no shared key to leak across the fleet.

## Authentication (GUI password)

The GUI has no login of its own; access control reuses the **kiosk password**.
**Implemented** in `server.py`: when a password-hash file exists and is
non-empty, every **non-localhost** request (LAN / Tailscale / remote) must log
in via `/login`; the walk-up kiosk on `127.0.0.1` is exempt -- physical presence
is the trust (override with `--require-local-auth`). No file -> no auth.

The hash file location is configurable -- `COMMEC_GUI_PASSWORD_FILE` (or
`--password-file`), **default `~/.config/commec-gui/password.hash`** so a
standalone (non-image) install needs no root. The box image can point it at a
system path (e.g. `/etc/commec-gui/password.hash`) if preferred. (A bad actor
who can rewrite that path/env is already inside the box, so the configurability
isn't a meaningful weakening.)

The first-boot setup that sets the kiosk user-account password must ALSO write
a **hash** of that same password to that file (you can't read the account
password back out of shadow, so capture it once and write both `passwd` and the
hash). Use the GUI env's werkzeug so the server can verify it:

    printf '%s' "$KIOSK_PASSWORD" | python -c \
      'import sys; from werkzeug.security import generate_password_hash; \
       print(generate_password_hash(sys.stdin.read()))' \
      > "$COMMEC_GUI_PASSWORD_FILE"    # mode 0600/0640, readable by the GUI user

Nothing is baked into the base image; the hash is written per-device at setup.
Auth is only meaningful over HTTPS (the server warns if a password is set but
TLS is off; the session cookie is marked Secure when serving HTTPS).

## Network hardening (inbound)

Clients want the kiosk on the internet for **outbound** fetches (e.g. pulling
sequences from Google Drive). Outbound is fine; the thing to control is
**inbound** reach to the GUI port.

- **Outbound: unrestricted** (kiosk browser -> Drive, etc.).
- **Inbound to the GUI port (`COMMEC_GUI_PORT`, default 8765): allow only**
  localhost, private LAN ranges, and Tailscale; DROP everything else (the
  public internet). NAT usually blocks inbound already, but don't rely on it --
  enforce a host firewall so a public IP or an accidental port-forward can't
  expose the GUI.
  - Allowed sources: `127.0.0.0/8`, `::1`, `10.0.0.0/8`, `172.16.0.0/12`,
    `192.168.0.0/16`, and **Tailscale** (`100.64.0.0/10`, or the `tailscale0`
    interface). Tailscale inbound is OK -- it's an authenticated overlay.
  - Example (nftables):
    ```
    table inet commec {
      chain input {
        type filter hook input priority 0; policy accept;
        iifname "tailscale0" tcp dport 8765 accept
        tcp dport 8765 ip saddr { 127.0.0.0/8, 10.0.0.0/8, 172.16.0.0/12, \
          192.168.0.0/16, 100.64.0.0/10 } accept
        tcp dport 8765 drop
      }
    }
    ```

This pairs with the GUI password: the firewall keeps the public internet from
reaching the port at all, and allowed LAN/Tailscale clients still authenticate.

## Third-party licensing

The GUI bundles one vendored asset: `gui/vendor/plotly-3.0.1.min.js`
(plotly.js, **MIT**, (c) Plotly Inc.) so the report renders offline. It ships
unmodified with its MIT notice in the file header -- see `gui/vendor/NOTICE.md`.
commec-gui is itself MIT, so this is compatible. **Carry this into the image's
and the upstream project's third-party / NOTICES accounting.**

## Storage & retention

Implemented entirely in `server.py` (`--runs-dir` / `--runs-keep` and the
periodic orphan-process sweep). There is **nothing for the image pipeline to do
here** -- no tmpfs mount, no swap-off. (This reverses an earlier
ephemeral/tmpfs design.)

Each run gets one persistent directory under `--runs-dir` (pointed at normal
persistent disk, e.g. somewhere under the kiosk user's home) holding the whole
run: `input.fasta`, `config.used.yaml`, commec's raw search intermediates
(`input_<prefix>/*.cleaned.fasta`/`.faa`/..., `output_<prefix>/*.blastn`/
`.blastx`/`.hmmscan`...), and the polished artifacts (`<prefix>.output.json`,
`<prefix>_summary.html`, `<prefix>.screen.log`). Nothing is wiped at end-of-run.

- **Capacity.** Disk now grows by roughly the full size of each run (bulky BLAST
  output included), so it accumulates rather than vanishing on power-off. Set
  `--runs-keep N` for a rolling cap (oldest whole run dirs pruned past the cap),
  or prune out of band. Sizing on the prototype is comfortable, but unbounded
  retention on a small disk will eventually fill it.
- **The sweep never deletes run dirs.** `--sweep-interval` reaps process groups
  orphaned by a crash/restart (a run that was live when the server died) and
  marks those runs **Interrupted** so the GUI can offer a one-click Resume
  (re-runs the cut-off step, re-uses the rest). Cleanup of finished runs is
  governed solely by `--runs-keep` and explicit deletes from the GUI.
