#!/bin/bash
# dry_run_test.sh -- Quick validation that all paths, configs, env vars,
# and Python imports resolve correctly WITHOUT running any benchmarks.
# Run this on the cluster BEFORE the full setup_and_run.sh to catch errors fast.
#
# Usage:  SEQPROC_DATA_DIR=/path/to/data \
#         SEQPROC_PROJECT_ROOT=/path/to/seqproc-paper-analysis \
#         bash scripts/dry_run_test.sh
#
# Or: source the env from setup_and_run.sh logic below.

set -euo pipefail

PASS=0
FAIL=0
WARN=0

pass() { echo "  [PASS] $1"; PASS=$((PASS+1)); }
fail() { echo "  [FAIL] $1"; FAIL=$((FAIL+1)); }
warn() { echo "  [WARN] $1"; WARN=$((WARN+1)); }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# ── Detect WORKDIR (same logic as setup_and_run.sh) ──
if [ -n "${WORKDIR:-}" ]; then
    : # user-provided
elif [ -d "/scratch0" ] && [ -w "/scratch0" ]; then
    WORKDIR="/scratch0/$USER/seqproc-bench"
elif [ -d "/scratch1" ] && [ -w "/scratch1" ]; then
    WORKDIR="/scratch1/$USER/seqproc-bench"
else
    WORKDIR="$HOME/seqproc-bench"
fi

# ── Activate micromamba bench env (same as setup_and_run.sh) ──
MICROMAMBA_ROOT="$HOME/.micromamba"
BENCH_BIN="$MICROMAMBA_ROOT/envs/bench/bin"
if [ -d "$BENCH_BIN" ]; then
    export PATH="$BENCH_BIN:$PATH"
fi

DATA_DIR="${SEQPROC_DATA_DIR:-$WORKDIR/data}"
PROJECT_ROOT="${SEQPROC_PROJECT_ROOT:-$ANALYSIS_ROOT}"
export SEQPROC_DATA_DIR="$DATA_DIR"
export SEQPROC_PROJECT_ROOT="$PROJECT_ROOT"

SEQPROC_BIN="${SEQPROC_BIN:-$WORKDIR/combine-lab/seqproc/target/release/seqproc}"
MATCHBOX_BIN="${MATCHBOX_BIN:-$WORKDIR/matchbox/target/release/matchbox}"
SPLITCODE_BIN="${SPLITCODE_BIN:-$WORKDIR/splitcode/build/src/splitcode}"
export SEQPROC_BIN MATCHBOX_BIN SPLITCODE_BIN

echo "============================================================"
echo "DRY RUN TEST -- Validating pipeline prerequisites"
echo "============================================================"
echo ""
echo "  WORKDIR:              $WORKDIR"
echo "  ANALYSIS_ROOT:        $ANALYSIS_ROOT"
echo "  SEQPROC_PROJECT_ROOT: $PROJECT_ROOT"
echo "  SEQPROC_DATA_DIR:     $DATA_DIR"
echo ""

# ── 1. Directories ──
echo "--- 1. Directories ---"
[ -d "$WORKDIR" ] && pass "WORKDIR exists: $WORKDIR" || fail "WORKDIR missing: $WORKDIR"
[ -d "$DATA_DIR" ] && pass "DATA_DIR exists: $DATA_DIR" || fail "DATA_DIR missing: $DATA_DIR"
[ -d "$PROJECT_ROOT" ] && pass "PROJECT_ROOT exists: $PROJECT_ROOT" || fail "PROJECT_ROOT missing: $PROJECT_ROOT"
[ -d "$PROJECT_ROOT/configs" ] && pass "configs/ dir exists" || fail "configs/ dir missing in $PROJECT_ROOT"
[ -d "$PROJECT_ROOT/scripts" ] && pass "scripts/ dir exists" || fail "scripts/ dir missing in $PROJECT_ROOT"
echo ""

# ── 2. Tool binaries ──
echo "--- 2. Tool binaries ---"
if [ -x "$SEQPROC_BIN" ]; then
    pass "seqproc binary: $SEQPROC_BIN"
    # Quick geom compile test (no data needed, just checks parsing)
    GEOM_FILE="$PROJECT_ROOT/configs/seqproc/splitseq_filter_edit.geom"
    if [ -f "$GEOM_FILE" ]; then
        # Run seqproc with --help to check it starts; actual geom test needs files
        if "$SEQPROC_BIN" --help >/dev/null 2>&1; then
            pass "seqproc binary executes"
        else
            fail "seqproc binary crashes on --help"
        fi
    fi
else
    fail "seqproc binary missing: $SEQPROC_BIN"
fi

[ -x "$MATCHBOX_BIN" ] && pass "matchbox binary: $MATCHBOX_BIN" || fail "matchbox binary missing: $MATCHBOX_BIN"
[ -x "$SPLITCODE_BIN" ] && pass "splitcode binary: $SPLITCODE_BIN" || fail "splitcode binary missing: $SPLITCODE_BIN"
echo ""

# ── 3. Config files referenced by DATASETS ──
echo "--- 3. Config files ---"
CONFIG_FILES=(
    "configs/seqproc/splitseq_filter_edit.geom"
    "configs/seqproc/splitseq_replacement_edit.geom"
    "configs/seqproc/splitseq_singleend_primer_edit.geom"
    "configs/seqproc/10x_v2.geom"
    "configs/seqproc/sciseq3_edit.geom"
    "configs/seqproc/splitseq_bc1_seq2seq.tsv"
    "configs/seqproc/splitseq_bc2_seq2seq.tsv"
    "configs/seqproc/splitseq_bc3_seq2seq.tsv"
    "configs/seqproc/splitseq_bc23_whitelist.txt"
    "configs/seqproc/splitseq_bc1_whitelist_6bp.txt"
    "configs/matchbox/splitseq_replacement.mb"
    "configs/matchbox/splitseq_singleend.mb"
    "configs/matchbox/10x_v2.mb"
    "configs/matchbox/sciseq3.mb"
    "configs/splitcode/splitseq_paper.config"
    "configs/splitcode/splitseq_singleend.config"
    "configs/splitcode/10x_v2.config"
    "configs/splitcode/sciseq3.config"
)
for cf in "${CONFIG_FILES[@]}"; do
    [ -f "$PROJECT_ROOT/$cf" ] && pass "$cf" || fail "MISSING: $PROJECT_ROOT/$cf"
done
echo ""

# ── 4. Geom files reference relative whitelist paths -- check them from PROJECT_ROOT ──
echo "--- 4. Whitelist paths in geom files (relative to PROJECT_ROOT) ---"
for geom in "$PROJECT_ROOT"/configs/seqproc/*.geom; do
    while IFS= read -r line; do
        # Extract quoted paths like "configs/seqproc/foo.txt"
        path=$(echo "$line" | grep -oP '"[^"]+\.(txt|tsv)"' | tr -d '"' || true)
        if [ -n "$path" ]; then
            if [ -f "$PROJECT_ROOT/$path" ]; then
                pass "$(basename "$geom"): $path"
            else
                fail "$(basename "$geom"): $path NOT FOUND (relative to $PROJECT_ROOT)"
            fi
        fi
    done < "$geom"
done
echo ""

# ── 5. Data files (full dataset) ──
echo "--- 5. Data files (full) ---"
DATA_FILES=(
    "SRR6750041_R1.fastq"
    "SRR6750041_R2.fastq"
    "SRR13948564_full.fastq"
    "10x_short/SRR8315379_R1.fastq"
    "10x_short/SRR8315379_R2.fastq"
    "SRR7827254_1.fastq"
    "SRR7827254_2.fastq"
)
for df in "${DATA_FILES[@]}"; do
    if [ -f "$DATA_DIR/$df" ]; then
        SIZE=$(du -h "$DATA_DIR/$df" | cut -f1)
        pass "$df ($SIZE)"
    else
        fail "MISSING: $DATA_DIR/$df"
    fi
done
echo ""

# ── 6. Python environment ──
echo "--- 6. Python environment ---"
if command -v python3 &>/dev/null; then
    PY_VER=$(python3 --version 2>&1)
    PY_PATH=$(which python3)
    pass "python3 found: $PY_VER ($PY_PATH)"
else
    fail "python3 not on PATH"
fi

# Check key imports
echo "  Checking Python imports..."
python3 -c "import matplotlib; print('    matplotlib', matplotlib.__version__)" 2>/dev/null && pass "matplotlib" || fail "matplotlib import"
python3 -c "import numpy; print('    numpy', numpy.__version__)" 2>/dev/null && pass "numpy" || fail "numpy import"
python3 -c "from dataclasses import dataclass; print('    dataclasses OK')" 2>/dev/null && pass "dataclasses" || fail "dataclasses import (Python too old?)"
echo ""

# ── 7. Python data_config resolution ──
echo "--- 7. data_config.py path resolution ---"
python3 -c "
import os, sys
sys.path.insert(0, '$PROJECT_ROOT/scripts')
os.environ['SEQPROC_DATA_DIR'] = '$DATA_DIR'
os.environ['SEQPROC_PROJECT_ROOT'] = '$PROJECT_ROOT'
from data_config import PROJECT_ROOT, DATA_DIR, CONFIGS, TOOL_CONFIGS
print(f'    PROJECT_ROOT = {PROJECT_ROOT}')
print(f'    DATA_DIR     = {DATA_DIR}')
print(f'    CONFIGS      = {CONFIGS}')
# Check a sample config file
sample = TOOL_CONFIGS['splitseq_pe']['seqproc_geom']
exists = sample.exists()
print(f'    Sample geom  = {sample}  exists={exists}')
if not exists:
    sys.exit(1)
" && pass "data_config paths resolve correctly" || fail "data_config path resolution broken"
echo ""

# ── 8. run_paper_benchmarks.py import + DATASETS check ──
echo "--- 8. run_paper_benchmarks.py DATASETS path check ---"
python3 -c "
import os, sys
sys.path.insert(0, '$PROJECT_ROOT/scripts')
os.environ['SEQPROC_DATA_DIR'] = '$DATA_DIR'
os.environ['SEQPROC_PROJECT_ROOT'] = '$PROJECT_ROOT'
os.environ['SEQPROC_BIN'] = '$SEQPROC_BIN'
os.environ['MATCHBOX_BIN'] = '$MATCHBOX_BIN'
os.environ['SPLITCODE_BIN'] = '$SPLITCODE_BIN'

# Import the module (checks all top-level code)
import run_paper_benchmarks as rpb

# Apply full reads level
rpb._apply_reads_level('full')

errors = 0
for key, ds in rpb.DATASETS.items():
    # Check config files
    for field in ['seqproc_geom', 'matchbox_config', 'splitcode_config', 'seqproc_geom_rev']:
        if field in ds and ds[field] is not None:
            p = ds[field]
            if not p.exists():
                print(f'    [FAIL] {key}.{field} = {p}')
                errors += 1
            else:
                print(f'    [OK]   {key}.{field}')
    # Check map files
    if 'seqproc_maps' in ds:
        for i, mp in enumerate(ds['seqproc_maps']):
            if not mp.exists():
                print(f'    [FAIL] {key}.seqproc_maps[{i}] = {mp}')
                errors += 1
            else:
                print(f'    [OK]   {key}.seqproc_maps[{i}]')
    # Check data files
    if ds.get('r1') and not ds['r1'].exists():
        print(f'    [FAIL] {key}.r1 = {ds[\"r1\"]}')
        errors += 1
    elif ds.get('r1'):
        print(f'    [OK]   {key}.r1')
    if ds.get('r2') and not ds['r2'].exists():
        print(f'    [FAIL] {key}.r2 = {ds[\"r2\"]}')
        errors += 1
    elif ds.get('r2'):
        print(f'    [OK]   {key}.r2')

if errors > 0:
    print(f'\\n    {errors} path(s) broken!')
    sys.exit(1)
else:
    print(f'\\n    All DATASETS paths OK')
" && pass "All DATASETS paths verified" || fail "DATASETS paths broken"
echo ""

# ── 9. run_all.sh sanity ──
echo "--- 9. run_all.sh ---"
[ -f "$PROJECT_ROOT/scripts/run_all.sh" ] && pass "run_all.sh exists" || fail "run_all.sh missing"
# Check it doesn't prefer venv python
if grep -q 'venv/bin/python3' "$PROJECT_ROOT/scripts/run_all.sh" 2>/dev/null; then
    fail "run_all.sh still references venv/bin/python3"
else
    pass "run_all.sh uses PATH python3 (no venv preference)"
fi
echo ""

# ── Summary ──
echo "============================================================"
echo "SUMMARY: $PASS passed, $FAIL failed, $WARN warnings"
echo "============================================================"
if [ "$FAIL" -gt 0 ]; then
    echo "FIX THE FAILURES ABOVE BEFORE RUNNING setup_and_run.sh"
    exit 1
else
    echo "All checks passed -- safe to run the full pipeline."
    exit 0
fi
