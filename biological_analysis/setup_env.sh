#!/usr/bin/env bash
# Create the locked downstream Python environment.
#   bash biological_analysis/setup_env.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
VENV=$HERE/.venv_downstream
UV=${UV_BIN:-uv}
PYTHON=${DOWNSTREAM_PYTHON:-/opt/local/stow/Python3-3.12.3/bin/python3}
CACHE=${DOWNSTREAM_UV_CACHE:-$HERE/work/uv-cache}
[ -x "$PYTHON" ] || { echo "Python 3.12 not found at $PYTHON" >&2; exit 2; }
command -v "$UV" >/dev/null || { echo "uv executable not found: $UV" >&2; exit 2; }
[ -f "$HERE/requirements.lock" ] || {
  echo "missing $HERE/requirements.lock; generate it from requirements.in with uv pip compile" >&2
  exit 2
}
mkdir -p "$CACHE"
UV_CACHE_DIR="$CACHE" "$UV" venv --clear --python "$PYTHON" "$VENV"
UV_CACHE_DIR="$CACHE" "$UV" pip sync --python "$VENV/bin/python" --require-hashes "$HERE/requirements.lock"
"$VENV/bin/python" -c "import scanpy, leidenalg, umap, sklearn, scipy, matplotlib; print('downstream env OK')"
echo "venv: $VENV"
