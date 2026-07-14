#!/usr/bin/env bash
# Generate a TLS certificate for serving commec-gui over the LAN, writing
# certs/server.crt and certs/server.key beside this script.
#
# Requires mkcert: the issued cert is trusted by browsers on machines where the
# mkcert root CA is installed (-> no warning). If mkcert is missing this errors
# out loudly (set COMMEC_TLS=0 to serve plain HTTP for local/dev instead) --
# we deliberately do NOT fall back to a self-signed openssl cert, so a missing
# mkcert is caught at provisioning rather than silently shipping warning-certs.
#
# The cert's Subject Alternative Names include every address clients use to
# reach the server (browsers ignore the legacy CN field): by default
# localhost, 127.0.0.1, this machine's hostname, and its detected LAN IP.
# Pass extra hostnames/IPs as arguments to add them.
#
#   ./gen-cert.sh                 # auto-detect
#   ./gen-cert.sh 192.168.1.42    # add an explicit address
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$HERE/certs"
mkdir -p "$OUT"
CRT="$OUT/server.crt"
KEY="$OUT/server.key"

HOSTNAME_S="$(hostname 2>/dev/null || echo localhost)"

# Best-effort LAN IPv4 detection (macOS first, then Linux).
LAN_IP=""
if command -v ipconfig >/dev/null 2>&1; then
  for IF in en0 en1 en2; do
    LAN_IP="$(ipconfig getifaddr "$IF" 2>/dev/null || true)"
    [[ -n "$LAN_IP" ]] && break
  done
fi
if [[ -z "$LAN_IP" ]] && command -v hostname >/dev/null 2>&1; then
  LAN_IP="$(hostname -I 2>/dev/null | awk '{print $1}' || true)"
fi

# Assemble the SAN list. NAMES is the flat list mkcert wants; DNS/IPS feed the
# openssl SAN string.
DNS=("localhost" "$HOSTNAME_S")
IPS=("127.0.0.1")
[[ -n "$LAN_IP" ]] && IPS+=("$LAN_IP")
for arg in "$@"; do
  if [[ "$arg" =~ ^[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+$ ]]; then IPS+=("$arg"); else DNS+=("$arg"); fi
done
NAMES=("${DNS[@]}" "${IPS[@]}")

if command -v mkcert >/dev/null 2>&1; then
  echo "Using mkcert for: ${NAMES[*]}"
  # Ensure the local CA is trusted on THIS machine. Adding it to the system
  # trust store needs sudo; if that fails (e.g. headless), the cert is still
  # issued -- browsers just show a warning until the CA is installed.
  if ! mkcert -install 2>/dev/null; then
    echo "NOTE: 'mkcert -install' could not update the system trust store"
    echo "      (needs an interactive sudo password). The cert is still valid;"
    echo "      run 'mkcert -install' yourself once to remove browser warnings."
  fi
  mkcert -cert-file "$CRT" -key-file "$KEY" "${NAMES[@]}"
  chmod 600 "$KEY"
  echo
  echo "mkcert root CA (install on LAN client machines to remove warnings):"
  echo "  $(mkcert -CAROOT)/rootCA.pem"
else
  echo "ERROR: mkcert not found." >&2
  echo "The GUI's TLS path requires mkcert. Install it (and run 'mkcert -install'" >&2
  echo "once so browsers trust the cert), or set COMMEC_TLS=0 in .env to serve" >&2
  echo "plain HTTP for local/dev use." >&2
  exit 1
fi

echo
echo "Wrote:"
echo "  $CRT"
echo "  $KEY"
