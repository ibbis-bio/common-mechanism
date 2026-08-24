#!/usr/bin/env bash
# Serve commec-gui to the LAN (bind 0.0.0.0). Reads config from an env file
# (default: .env beside this script; override by passing a path as $1).
#
# If COMMEC_TLS_CERT and COMMEC_TLS_KEY are set, serves over HTTPS. Generate a
# self-signed pair with ./gen-cert.sh first if you don't have one.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${1:-$HERE/.env}"

if [[ -f "$ENV_FILE" ]]; then
  set -a
  # shellcheck disable=SC1090
  source "$ENV_FILE"
  set +a
else
  echo "No env file at $ENV_FILE (copy .env.example to .env)." >&2
fi

COMMEC_CONDA_ENV="${COMMEC_CONDA_ENV:-commec-dev}"
COMMEC_GUI_PORT="${COMMEC_GUI_PORT:-443}"

if [[ "${CONDA_DEFAULT_ENV:-}" != "$COMMEC_CONDA_ENV" ]]; then
  if ! CONDA_BASE="$(conda info --base 2>/dev/null)"; then
    echo "conda not found on PATH; cannot activate '$COMMEC_CONDA_ENV'." >&2
    exit 1
  fi
  # shellcheck disable=SC1091
  source "$CONDA_BASE/etc/profile.d/conda.sh"
  conda activate "$COMMEC_CONDA_ENV"
fi

ARGS=(--lan --port "$COMMEC_GUI_PORT")
if [[ -n "${COMMEC_DB_DIR:-}" ]]; then
  ARGS+=(--databases "$COMMEC_DB_DIR")
fi
if [[ -n "${COMMEC_GUI_THREADS:-}" ]]; then
  ARGS+=(--threads "$COMMEC_GUI_THREADS")
fi
if [[ -n "${COMMEC_GUI_PASSWORD_FILE:-}" ]]; then
  ARGS+=(--password-file "$COMMEC_GUI_PASSWORD_FILE")
fi

# TLS: use an explicit cert/key pair if configured, otherwise --tls-auto
# generates one on first startup. Set COMMEC_TLS=0 to force plain HTTP.
if [[ "${COMMEC_TLS:-1}" == "0" ]]; then
  echo "WARNING: COMMEC_TLS=0 -- serving LAN traffic over plain HTTP." >&2
elif [[ -n "${COMMEC_TLS_CERT:-}" || -n "${COMMEC_TLS_KEY:-}" ]]; then
  if [[ -z "${COMMEC_TLS_CERT:-}" || -z "${COMMEC_TLS_KEY:-}" ]]; then
    echo "Set BOTH COMMEC_TLS_CERT and COMMEC_TLS_KEY, or neither." >&2
    exit 1
  fi
  ARGS+=(--tls-cert "$COMMEC_TLS_CERT" --tls-key "$COMMEC_TLS_KEY")
else
  ARGS+=(--tls-auto)
fi

exec commec gui "${ARGS[@]}"
