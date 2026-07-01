#!/usr/bin/env bash
# One-time setup of the Phase 2A python env (scanpy + clustering + notebook tooling).
# Self-contained venv so it reproduces identically on the cluster regardless of system packages.
#   bash biological_analysis/setup_env.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
VENV=$HERE/.venv_phase2a
# --copies: module-provided pythons on HPC build broken symlinks into a venv, so copy the
# interpreter in. Drive pip via "python -m pip" (robust even if the bin/pip shebang is off).
rm -rf "$VENV"
python3 -m venv --copies "$VENV"
"$VENV/bin/python" -m pip install --quiet --upgrade pip wheel
"$VENV/bin/python" -m pip install --quiet \
  "scanpy==1.12.2" "anndata>=0.10" leidenalg igraph "scikit-learn>=1.3" \
  "scipy>=1.11" "numpy>=1.26,<3" "matplotlib>=3.6" "pandas>=2.0" \
  nbformat ipykernel nbclient
"$VENV/bin/python" -c "import scanpy, leidenalg, umap, sklearn, scipy, matplotlib; print('phase2a env OK')"
echo "venv: $VENV"
