#!/usr/bin/env bash
# Spin up commec-gui in kiosk mode: bind localhost and auto-open a browser.
# Configuration is read from an env file (default: .env beside this script;
# override by passing a path as the first argument).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${1:-$HERE/.env}"

if [[ -f "$ENV_FILE" ]]; then
  set -a            # export everything defined while sourcing
  # shellcheck disable=SC1090
  source "$ENV_FILE"
  set +a
else
  echo "No env file at $ENV_FILE (copy .env.example to .env)." >&2
  echo "Continuing with defaults; no database dir will be prefilled." >&2
fi

COMMEC_CONDA_ENV="${COMMEC_CONDA_ENV:-commec-dev}"
COMMEC_GUI_PORT="${COMMEC_GUI_PORT:-8765}"

# Activate the conda env that has commec + flask, unless already active.
if [[ "${CONDA_DEFAULT_ENV:-}" != "$COMMEC_CONDA_ENV" ]]; then
  if ! CONDA_BASE="$(conda info --base 2>/dev/null)"; then
    echo "conda not found on PATH; cannot activate '$COMMEC_CONDA_ENV'." >&2
    exit 1
  fi
  # shellcheck disable=SC1091
  source "$CONDA_BASE/etc/profile.d/conda.sh"
  conda activate "$COMMEC_CONDA_ENV"
fi

ARGS=(--kiosk --port "$COMMEC_GUI_PORT")
if [[ -n "${COMMEC_DB_DIR:-}" ]]; then
  ARGS+=(--databases "$COMMEC_DB_DIR")
fi
if [[ -n "${COMMEC_GUI_THREADS:-}" ]]; then
  ARGS+=(--threads "$COMMEC_GUI_THREADS")
fi

exec python "$HERE/server.py" "${ARGS[@]}"
