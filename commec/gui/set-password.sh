#!/usr/bin/env bash
# Create or replace the Commec GUI access-password hash. This is intentionally
# a GUI-local helper, not a `commec` subcommand: it handles a deployment secret,
# not screening. It reads the same .env file as kiosk.sh and lan.sh.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${1:-$HERE/.env}"

if [[ -f "$ENV_FILE" ]]; then
  set -a
  # shellcheck disable=SC1090
  source "$ENV_FILE"
  set +a
elif [[ $# -gt 0 ]]; then
  echo "Environment file not found: $ENV_FILE" >&2
  exit 1
fi

COMMEC_CONDA_ENV="${COMMEC_CONDA_ENV:-commec-dev}"
HASH_FILE="${COMMEC_GUI_PASSWORD_FILE:-$HOME/.config/commec-gui/password.hash}"

# Match the launchers: use the GUI environment's Werkzeug, rather than relying
# on whatever Python happens to be on the caller's PATH.
if [[ "${CONDA_DEFAULT_ENV:-}" != "$COMMEC_CONDA_ENV" ]]; then
  if ! CONDA_BASE="$(conda info --base)"; then
    echo "conda not found on PATH; cannot activate '$COMMEC_CONDA_ENV'." >&2
    exit 1
  fi
  # shellcheck disable=SC1091
  source "$CONDA_BASE/etc/profile.d/conda.sh"
  conda activate "$COMMEC_CONDA_ENV"
fi

read -r -s -p "New GUI password: " password
printf '\n' >&2
read -r -s -p "Confirm GUI password: " confirm
printf '\n' >&2

if [[ -z "$password" ]]; then
  echo "Password cannot be empty." >&2
  exit 1
fi
if [[ ${#password} -lt 8 ]]; then
  echo "Use at least 8 characters." >&2
  exit 1
fi
if [[ "$password" != "$confirm" ]]; then
  echo "Passwords did not match." >&2
  exit 1
fi
unset confirm

hash_dir="$(dirname "$HASH_FILE")"
install -d -m 700 "$hash_dir"
umask 077
tmp="$(mktemp "$hash_dir/.password.hash.XXXXXX")"
trap 'rm -f "$tmp"' EXIT

printf '%s' "$password" | python -c \
  'import sys; from werkzeug.security import generate_password_hash; print(generate_password_hash(sys.stdin.read()))' \
  > "$tmp"
unset password
chmod 600 "$tmp"
mv -f "$tmp" "$HASH_FILE"
trap - EXIT

echo "GUI password hash written to $HASH_FILE. Restart the GUI to apply it."
