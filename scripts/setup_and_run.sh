#!/bin/bash
set -euo pipefail

###############################################################################
# seqproc Full Paper Benchmark -- Clean Machine Setup + Run
#
# NO SUDO REQUIRED. Installs everything to user-local directories.
# Assumes: Linux x86_64, internet access, git, curl, gcc/g++, make.
# Python and cmake are installed automatically if missing or too old.
#
# Every step is guarded: re-running this script skips completed work.
#
# Usage (from the cloned analysis repo root):
#   chmod +x scripts/setup_and_run.sh && ./scripts/setup_and_run.sh
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# WORKDIR: override with env var, or auto-detect project/scratch space, or fall back to $HOME
if [ -n "${WORKDIR:-}" ]; then
    : # user-provided, keep it
elif [ -d "/scratch0" ] && [ -w "/scratch0" ]; then
    WORKDIR="/scratch0/$USER/seqproc-bench"
elif [ -d "/scratch1" ] && [ -w "/scratch1" ]; then
    WORKDIR="/scratch1/$USER/seqproc-bench"
else
    WORKDIR="$HOME/seqproc-bench"
fi
LOCAL_BIN="$HOME/.local/bin"
MICROMAMBA_ROOT="$HOME/.micromamba"
THREADS=$(nproc 2>/dev/null || echo 4)
REPLICATES=3

SEQPROC_BIN="$WORKDIR/combine-lab/seqproc/target/release/seqproc"
MATCHBOX_BIN="$WORKDIR/matchbox/target/release/matchbox"
SPLITCODE_BIN="$WORKDIR/splitcode/build/src/splitcode"

mkdir -p "$WORKDIR" "$LOCAL_BIN"
export PATH="$LOCAL_BIN:$HOME/.cargo/bin:$PATH"

# Tee all output to a timestamped log file AND the terminal
LOG_FILE="$WORKDIR/setup_and_run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "Logging to: $LOG_FILE"

# GitHub auth: use GH_TOKEN for HTTPS cloning of private repos.
# If not set, fall back to SSH (requires SSH keys).
if [ -n "${GH_TOKEN:-}" ]; then
    GH_PREFIX="https://${GH_TOKEN}@github.com/"
else
    GH_PREFIX="git@github.com:"
fi

echo "================================================================"
echo "seqproc Full Paper Benchmark Setup"
echo "  WORKDIR:       $WORKDIR"
echo "  ANALYSIS_ROOT: $ANALYSIS_ROOT"
echo "  THREADS:       $THREADS"
echo "  REPLICATES:    $REPLICATES"
echo "================================================================"

# ---------------------------------------------------------------------------
# Helper: put micromamba bench env bin/ on PATH (no activate needed)
# ---------------------------------------------------------------------------
BENCH_BIN="$MICROMAMBA_ROOT/envs/bench/bin"
activate_bench_env() {
    if [ -d "$BENCH_BIN" ]; then
        export PATH="$BENCH_BIN:$LOCAL_BIN:$HOME/.cargo/bin:$PATH"
    else
        export PATH="$LOCAL_BIN:$HOME/.cargo/bin:$PATH"
    fi
}

# ---------------------------------------------------------------------------
# Helper: ensure micromamba is installed
# ---------------------------------------------------------------------------
ensure_micromamba() {
    if [ ! -x "$MICROMAMBA_ROOT/bin/micromamba" ]; then
        echo "  Installing micromamba..."
        mkdir -p "$MICROMAMBA_ROOT/bin"
        MAMBA_DL="$WORKDIR/_micromamba_download.tar.bz2"
        curl -fSL -o "$MAMBA_DL" https://micro.mamba.pm/api/micromamba/linux-64/latest
        if ! file "$MAMBA_DL" | grep -q "bzip2"; then
            echo "  [ERROR] micromamba download is not a valid bzip2 file."
            echo "  Contents: $(head -c 200 "$MAMBA_DL")"
            echo "  Try downloading manually from https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html"
            rm -f "$MAMBA_DL"
            exit 1
        fi
        tar -xjf "$MAMBA_DL" -C "$MICROMAMBA_ROOT/bin" --strip-components=1 bin/micromamba
        chmod +x "$MICROMAMBA_ROOT/bin/micromamba"
        rm -f "$MAMBA_DL"
    fi
    export MAMBA_ROOT_PREFIX="$MICROMAMBA_ROOT"
}

# ============================================================================
# 1. Check / install prerequisites (no sudo)
# ============================================================================
echo "[1/8] Checking prerequisites..."

require_cmd() {
    if ! command -v "$1" &>/dev/null; then
        echo "  [ERROR] Required command '$1' not found."
        echo "  On an HPC cluster, try: module load $1"
        echo "  Or ask your sysadmin to install: $2"
        exit 1
    fi
    echo "  [OK] $1"
}

require_cmd git      "git"
require_cmd curl     "curl"
require_cmd gcc      "gcc / build-essential"
require_cmd make     "make / build-essential"

# -- Python >= 3.9 (install via micromamba if missing/old) --
NEED_MAMBA_PYTHON=false
if command -v python3 &>/dev/null; then
    PY_VER=$(python3 -c 'import sys; print(sys.version_info.minor)' 2>/dev/null || echo 0)
    if [ "$PY_VER" -ge 9 ] 2>/dev/null; then
        echo "  [OK] python3 ($(python3 --version))"
    else
        NEED_MAMBA_PYTHON=true
    fi
else
    NEED_MAMBA_PYTHON=true
fi

if [ "$NEED_MAMBA_PYTHON" = true ]; then
    echo "  Python >= 3.9 not found."
    ensure_micromamba
    if [ ! -d "$MICROMAMBA_ROOT/envs/bench" ]; then
        echo "  Creating micromamba bench env (python 3.11 + sra-tools)..."
        "$MICROMAMBA_ROOT/bin/micromamba" create -y -n bench \
            python=3.11 "numpy=1.26.*" "matplotlib=3.8.*" sra-tools \
            -c conda-forge -c bioconda
    fi
    activate_bench_env
    echo "  [OK] python3 via micromamba ($($BENCH_BIN/python3 --version))"
fi

# -- cmake (install locally if missing) --
if ! command -v cmake &>/dev/null; then
    if [ -x "$LOCAL_BIN/cmake" ]; then
        echo "  [OK] cmake (local)"
    else
        echo "  Installing cmake locally..."
        CMAKE_VER="3.28.3"
        curl -sL "https://github.com/Kitware/CMake/releases/download/v${CMAKE_VER}/cmake-${CMAKE_VER}-linux-x86_64.tar.gz" \
            | tar -xz -C "$WORKDIR"
        ln -sf "$WORKDIR/cmake-${CMAKE_VER}-linux-x86_64/bin/cmake" "$LOCAL_BIN/cmake"
        echo "  [OK] cmake (local install)"
    fi
else
    echo "  [OK] cmake"
fi

# -- Rust toolchain --
if ! command -v cargo &>/dev/null; then
    echo "  Installing Rust toolchain..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --no-modify-path
    source "$HOME/.cargo/env"
fi
export PATH="$HOME/.cargo/bin:$PATH"
echo "  [OK] cargo $(cargo --version 2>/dev/null | head -1)"

# -- SRA toolkit (via micromamba if not already available) --
activate_bench_env
if ! command -v fasterq-dump &>/dev/null; then
    echo "  fasterq-dump not found. Installing sra-tools via micromamba..."
    ensure_micromamba
    if [ ! -d "$MICROMAMBA_ROOT/envs/bench" ]; then
        "$MICROMAMBA_ROOT/bin/micromamba" create -y -n bench sra-tools -c conda-forge -c bioconda
    else
        "$MICROMAMBA_ROOT/bin/micromamba" install -y -n bench sra-tools -c conda-forge -c bioconda
    fi
    activate_bench_env
    echo "  [OK] fasterq-dump via micromamba"
else
    echo "  [OK] fasterq-dump"
fi

echo ""

# ============================================================================
# 2. Clone and build seqproc + ANTISEQUENCE
# ============================================================================
echo "[2/8] seqproc + ANTISEQUENCE..."
mkdir -p "$WORKDIR/combine-lab"

if [ ! -d "$WORKDIR/combine-lab/ANTISEQUENCE" ]; then
    git clone --branch cleanup_and_final_touches \
        "${GH_PREFIX}COMBINE-lab/ANTISEQUENCE.git" \
        "$WORKDIR/combine-lab/ANTISEQUENCE"
else
    echo "  Updating ANTISEQUENCE..."
    git -C "$WORKDIR/combine-lab/ANTISEQUENCE" pull --ff-only 2>/dev/null || true
fi

if [ ! -d "$WORKDIR/combine-lab/seqproc" ]; then
    git clone --branch edit_distance_map \
        "${GH_PREFIX}COMBINE-lab/seqproc.git" \
        "$WORKDIR/combine-lab/seqproc"
else
    echo "  Updating seqproc..."
    git -C "$WORKDIR/combine-lab/seqproc" pull --ff-only 2>/dev/null || true
fi

SEQPROC_STAMP="$WORKDIR/combine-lab/seqproc/.build_commit"
SEQPROC_HEAD="$(git -C "$WORKDIR/combine-lab/seqproc" rev-parse HEAD 2>/dev/null || echo none)"
SEQPROC_BUILT="$(cat "$SEQPROC_STAMP" 2>/dev/null || echo missing)"
if [ -x "$SEQPROC_BIN" ] && [ "$SEQPROC_HEAD" = "$SEQPROC_BUILT" ]; then
    echo "  [SKIP] seqproc already built at $SEQPROC_HEAD"
else
    echo "  Building seqproc (release) at $SEQPROC_HEAD..."
    cd "$WORKDIR/combine-lab/seqproc"
    cargo build --release
    echo "$SEQPROC_HEAD" > "$SEQPROC_STAMP"
fi
echo "  seqproc: $("$SEQPROC_BIN" --version 2>&1 || echo 'built')"

# ============================================================================
# 3. Clone and build matchbox
# ============================================================================
echo "[3/8] matchbox..."
if [ ! -d "$WORKDIR/matchbox" ]; then
    git clone https://github.com/jakob-schuster/matchbox.git \
        "$WORKDIR/matchbox"
fi

if [ -x "$MATCHBOX_BIN" ]; then
    echo "  [SKIP] matchbox already built: $MATCHBOX_BIN"
else
    echo "  Building matchbox (release)..."
    cd "$WORKDIR/matchbox"
    cargo build --release
fi

# ============================================================================
# 4. Clone and build splitcode
# ============================================================================
echo "[4/8] splitcode..."
if [ ! -d "$WORKDIR/splitcode" ]; then
    git clone https://github.com/pachterlab/splitcode.git \
        "$WORKDIR/splitcode"
fi

if [ -x "$SPLITCODE_BIN" ]; then
    echo "  [SKIP] splitcode already built: $SPLITCODE_BIN"
else
    echo "  Building splitcode..."
    cd "$WORKDIR/splitcode"
    # GCC < 9 needs -lstdc++fs for std::filesystem
    GCC_MAJOR=$(gcc -dumpversion | cut -d. -f1)
    if [ "$GCC_MAJOR" -lt 9 ] 2>/dev/null; then
        if ! grep -q "stdc++fs" src/CMakeLists.txt; then
            echo "  Patching splitcode for GCC $GCC_MAJOR (adding -lstdc++fs)..."
            sed -i 's/target_link_libraries(splitcode splitcode_core pthread)/target_link_libraries(splitcode splitcode_core pthread stdc++fs)/' src/CMakeLists.txt
        fi
    fi
    rm -rf build && mkdir -p build && cd build
    cmake .. && make -j"$THREADS"
fi

# ============================================================================
# 5. Set up analysis repo + Python deps
# ============================================================================
echo "[5/8] Analysis repo + Python deps..."

if [ -f "$ANALYSIS_ROOT/scripts/run_all.sh" ]; then
    echo "  Running from cloned repo at $ANALYSIS_ROOT"
    if [ "$ANALYSIS_ROOT" != "$WORKDIR/seqproc-paper-analysis" ]; then
        ln -sfn "$ANALYSIS_ROOT" "$WORKDIR/seqproc-paper-analysis"
    fi
else
    if [ ! -d "$WORKDIR/seqproc-paper-analysis" ]; then
        git clone --branch phase3-orientation-benchmarks \
            "${GH_PREFIX}COMBINE-lab/seqproc-paper-analysis.git" \
            "$WORKDIR/seqproc-paper-analysis"
    fi
fi

cd "$WORKDIR/seqproc-paper-analysis"
activate_bench_env

# Install Python deps into micromamba bench env
# numpy/matplotlib MUST come from conda-forge with numpy<2.0 (numpy 2.x
# requires X86_V2 CPU instructions that older cluster nodes lack).
if [ -d "$BENCH_BIN" ]; then
    # Only fix numpy/matplotlib if they fail to import (avoid breaking a working env)
    if ! "$BENCH_BIN/python3" -c "import numpy, matplotlib" 2>/dev/null; then
        echo "  Fixing numpy/matplotlib (need conda-forge numpy<2.0)..."
        "$BENCH_BIN/pip" uninstall -y numpy matplotlib 2>/dev/null || true
        rm -rf "$MICROMAMBA_ROOT/envs/bench/lib/python3.11/site-packages/numpy"* 2>/dev/null || true
        "$MICROMAMBA_ROOT/bin/micromamba" install -y -n bench "numpy=1.26.*" "matplotlib=3.8.*" \
            -c conda-forge --force-reinstall
    else
        echo "  [OK] numpy/matplotlib already working"
    fi
    # Install remaining pip deps (requirements.txt has no numpy/matplotlib)
    "$BENCH_BIN/pip" install -q -r requirements.txt 2>/dev/null \
        || "$BENCH_BIN/pip" install -r requirements.txt
    echo "  [OK] Python deps installed in micromamba bench env"
else
    echo "  [ERROR] micromamba bench env not found at $BENCH_BIN"
    echo "  Cannot install Python deps. Ensure step 1 completed."
    exit 1
fi

# ============================================================================
# 6. Download FULL SRA datasets
# ============================================================================
echo "[6/8] Downloading full SRA datasets..."
activate_bench_env
SRA_TMP="$WORKDIR/tmp"
DATA_DIR="$WORKDIR/data"
mkdir -p "$DATA_DIR/10x_short" "$SRA_TMP"

# Redirect SRA cache to project space (avoids NFS quota)
export NCBI_SETTINGS="$WORKDIR/ncbi_settings.kfg"
if [ ! -f "$NCBI_SETTINGS" ]; then
    mkdir -p "$WORKDIR/ncbi_cache"
    printf '/repository/user/main/public/root = "%s"\n' "$WORKDIR/ncbi_cache" > "$NCBI_SETTINGS"
fi

# Tell Python scripts where data and configs live (data_config.py reads these)
export SEQPROC_DATA_DIR="$DATA_DIR"
export SEQPROC_PROJECT_ROOT="$ANALYSIS_ROOT"
echo "  SEQPROC_DATA_DIR=$DATA_DIR"
echo "  SEQPROC_PROJECT_ROOT=$ANALYSIS_ROOT"
cd "$DATA_DIR"

if [ ! -f SRR6750041_R1.fastq ]; then
    echo "  Downloading SRR6750041 (SPLiT-seq PE, ~20 GB)..."
    fasterq-dump --split-files SRR6750041 --threads "$THREADS" --temp "$SRA_TMP" --outdir "$DATA_DIR"
    mv SRR6750041_1.fastq SRR6750041_R1.fastq
    mv SRR6750041_2.fastq SRR6750041_R2.fastq
    rm -rf SRR6750041/
    head -40000000 SRR6750041_R1.fastq > SRR6750041_10M_R1.fastq
    head -40000000 SRR6750041_R2.fastq > SRR6750041_10M_R2.fastq
else
    echo "  [SKIP] SRR6750041 already present"
fi

if [ ! -f SRR13948564_full.fastq ]; then
    echo "  Downloading SRR13948564 (LR-SPLiT-seq, ~5 GB)..."
    fasterq-dump SRR13948564 --threads "$THREADS" --temp "$SRA_TMP" --outdir "$DATA_DIR"
    mv SRR13948564.fastq SRR13948564_full.fastq
    rm -rf SRR13948564/
else
    echo "  [SKIP] SRR13948564 already present"
fi

if [ ! -f 10x_short/SRR8315379_R1.fastq ]; then
    echo "  Downloading SRR8315379 (10x Chromium v2, ~10 GB)..."
    fasterq-dump --split-files SRR8315379 --threads "$THREADS" --temp "$SRA_TMP" --outdir "$DATA_DIR"
    mv SRR8315379_1.fastq 10x_short/SRR8315379_R1.fastq
    mv SRR8315379_2.fastq 10x_short/SRR8315379_R2.fastq
    rm -rf SRR8315379/
else
    echo "  [SKIP] SRR8315379 already present"
fi

if [ ! -f SRR7827254_1.fastq ]; then
    echo "  Downloading SRR7827254 (sci-RNA-seq3, ~3 GB)..."
    fasterq-dump --split-files SRR7827254 --threads "$THREADS" --temp "$SRA_TMP" --outdir "$DATA_DIR"
    rm -rf SRR7827254/
else
    echo "  [SKIP] SRR7827254 already present"
fi

cd "$WORKDIR/seqproc-paper-analysis"
echo "  Verifying data availability..."
python scripts/data_config.py --reads full

# ============================================================================
# 7. Run the full benchmark pipeline
# ============================================================================
echo "[7/8] Running full benchmark pipeline..."
activate_bench_env
export SEQPROC_BIN MATCHBOX_BIN SPLITCODE_BIN

chmod +x scripts/run_all.sh
bash scripts/run_all.sh --reads full --threads "$THREADS" --replicates "$REPLICATES"

# ============================================================================
# 8. Summary
# ============================================================================
echo ""
echo "================================================================"
echo "[8/8] DONE"
echo "================================================================"
echo ""
echo "Results:  $WORKDIR/seqproc-paper-analysis/results/paper_figures/"
echo ""
echo "Figures generated:"
ls -lh "$WORKDIR/seqproc-paper-analysis/results/paper_figures/"*.pdf 2>/dev/null || echo "  (check results directory)"
echo ""
echo "JSON data:"
ls -lh "$WORKDIR/seqproc-paper-analysis/results/paper_figures/"*.json 2>/dev/null || echo "  (check results directory)"
echo ""
echo "To copy figures into the paper repo:"
echo "  cp $WORKDIR/seqproc-paper-analysis/results/paper_figures/*.pdf /path/to/paper/Figures/"
