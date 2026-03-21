#!/bin/bash
set -euo pipefail

###############################################################################
# seqproc Full Paper Benchmark -- Clean Machine Setup + Run
#
# NO SUDO REQUIRED. Installs everything to $HOME/.local and $HOME/.cargo.
# Assumes: Linux x86_64, internet access, git, curl, python3, gcc/g++, cmake.
#
# Usage (from the cloned analysis repo root):
#   chmod +x scripts/setup_and_run.sh && ./scripts/setup_and_run.sh
#
# Cluster usage (typical HPC with module system):
#   module load gcc cmake python3 git curl   # load whatever is available
#   ./scripts/setup_and_run.sh
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
WORKDIR="$HOME/seqproc-bench"
LOCAL_BIN="$HOME/.local/bin"
THREADS=$(nproc 2>/dev/null || echo 4)
REPLICATES=3

mkdir -p "$WORKDIR" "$LOCAL_BIN"
export PATH="$LOCAL_BIN:$HOME/.cargo/bin:$PATH"

echo "================================================================"
echo "seqproc Full Paper Benchmark Setup"
echo "  WORKDIR:       $WORKDIR"
echo "  ANALYSIS_ROOT: $ANALYSIS_ROOT"
echo "  LOCAL_BIN:     $LOCAL_BIN"
echo "  THREADS:       $THREADS"
echo "  REPLICATES:    $REPLICATES"
echo "================================================================"

# ============================================================================
# 1. Check / install prerequisites (no sudo)
# ============================================================================
echo "[1/8] Checking prerequisites..."

# Helper: bail with a message if a required system tool is missing
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

# Python >= 3.9 is required for matplotlib >= 3.7 and numpy >= 1.24.
# If the system python3 is too old (or missing), install Python 3.11
# via micromamba (no sudo needed).
MICROMAMBA_ROOT="$WORKDIR/micromamba"
PYTHON_BIN=""
if command -v python3 &>/dev/null; then
    PY_VER=$(python3 -c 'import sys; print(sys.version_info.minor)' 2>/dev/null || echo 0)
    if [ "$PY_VER" -ge 9 ] 2>/dev/null; then
        PYTHON_BIN="python3"
        echo "  [OK] python3 ($(python3 --version))"
    fi
fi
if [ -z "$PYTHON_BIN" ]; then
    echo "  Python >= 3.9 not found. Installing Python 3.11 via micromamba..."
    if [ ! -x "$MICROMAMBA_ROOT/bin/micromamba" ]; then
        mkdir -p "$MICROMAMBA_ROOT/bin"
        curl -sL https://micro.mamba.pm/api/micromamba/linux-64/latest \
            | tar -xj -C "$MICROMAMBA_ROOT/bin" --strip-components=1 bin/micromamba
        chmod +x "$MICROMAMBA_ROOT/bin/micromamba"
    fi
    export MAMBA_ROOT_PREFIX="$MICROMAMBA_ROOT/envs"
    if [ ! -d "$MAMBA_ROOT_PREFIX/envs/bench" ]; then
        "$MICROMAMBA_ROOT/bin/micromamba" create -y -n bench python=3.11 -c conda-forge
    fi
    eval "$("$MICROMAMBA_ROOT/bin/micromamba" shell hook -s bash)"
    micromamba activate bench
    PYTHON_BIN="python3"
    echo "  [OK] python3 via micromamba ($(python3 --version))"
fi

# cmake -- try to find it, or install locally
if ! command -v cmake &>/dev/null; then
    echo "  cmake not found, installing locally..."
    CMAKE_VER="3.28.3"
    curl -sL "https://github.com/Kitware/CMake/releases/download/v${CMAKE_VER}/cmake-${CMAKE_VER}-linux-x86_64.tar.gz" \
        | tar -xz -C "$WORKDIR"
    ln -sf "$WORKDIR/cmake-${CMAKE_VER}-linux-x86_64/bin/cmake" "$LOCAL_BIN/cmake"
    echo "  [OK] cmake (local install)"
else
    echo "  [OK] cmake"
fi

# Rust toolchain (user-local, no sudo)
if ! command -v cargo &>/dev/null; then
    echo "  Installing Rust toolchain (user-local)..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --no-modify-path
    source "$HOME/.cargo/env"
fi
export PATH="$HOME/.cargo/bin:$PATH"
echo "  [OK] cargo $(cargo --version 2>/dev/null | head -1)"

# SRA toolkit (user-local, no sudo)
if ! command -v fasterq-dump &>/dev/null; then
    echo "  Installing SRA toolkit (user-local)..."
    SRA_TMP="$WORKDIR/_sra_install"
    mkdir -p "$SRA_TMP"
    curl -L -o "$SRA_TMP/sratoolkit.tar.gz" \
        https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
    tar -xzf "$SRA_TMP/sratoolkit.tar.gz" -C "$SRA_TMP"
    SRA_DIR=$(ls -d "$SRA_TMP"/sratoolkit.* 2>/dev/null | head -1)
    cp "$SRA_DIR/bin/fasterq-dump" "$LOCAL_BIN/"
    cp "$SRA_DIR/bin/prefetch" "$LOCAL_BIN/"
    chmod +x "$LOCAL_BIN/fasterq-dump" "$LOCAL_BIN/prefetch"
    rm -rf "$SRA_TMP"
    echo "  [OK] fasterq-dump (local install)"
else
    echo "  [OK] fasterq-dump"
fi

echo ""

# ============================================================================
# 2. Clone and build seqproc + ANTISEQUENCE
# ============================================================================
echo "[2/8] Cloning and building seqproc + ANTISEQUENCE..."
mkdir -p "$WORKDIR/combine-lab"

# ANTISEQUENCE (must be cloned first -- seqproc depends on it via path)
if [ ! -d "$WORKDIR/combine-lab/ANTISEQUENCE" ]; then
    git clone --branch cleanup_and_final_touches \
        git@github.com:COMBINE-lab/ANTISEQUENCE.git \
        "$WORKDIR/combine-lab/ANTISEQUENCE"
fi

# seqproc
if [ ! -d "$WORKDIR/combine-lab/seqproc" ]; then
    git clone --branch edit_distance_map \
        git@github.com:COMBINE-lab/seqproc.git \
        "$WORKDIR/combine-lab/seqproc"
fi

# Build seqproc (release mode)
cd "$WORKDIR/combine-lab/seqproc"
cargo build --release
SEQPROC_BIN="$WORKDIR/combine-lab/seqproc/target/release/seqproc"
echo "  seqproc binary: $SEQPROC_BIN"
"$SEQPROC_BIN" --version || true

# ============================================================================
# 3. Clone and build matchbox
# ============================================================================
echo "[3/8] Cloning and building matchbox..."
if [ ! -d "$WORKDIR/matchbox" ]; then
    git clone https://github.com/jakob-schuster/matchbox.git \
        "$WORKDIR/matchbox"
fi
cd "$WORKDIR/matchbox"
cargo build --release
MATCHBOX_BIN="$WORKDIR/matchbox/target/release/matchbox"
echo "  matchbox binary: $MATCHBOX_BIN"

# ============================================================================
# 4. Clone and build splitcode
# ============================================================================
echo "[4/8] Cloning and building splitcode..."
if [ ! -d "$WORKDIR/splitcode" ]; then
    git clone https://github.com/pachterlab/splitcode.git \
        "$WORKDIR/splitcode"
fi
cd "$WORKDIR/splitcode"
# GCC < 9 needs -lstdc++fs for std::filesystem support.
# CMAKE_EXE_LINKER_FLAGS does not propagate to the splitcode target,
# so we patch target_link_libraries in src/CMakeLists.txt directly.
GCC_MAJOR=$(gcc -dumpversion | cut -d. -f1)
if [ "$GCC_MAJOR" -lt 9 ] 2>/dev/null; then
    if ! grep -q "stdc++fs" src/CMakeLists.txt; then
        echo "  Patching splitcode for GCC $GCC_MAJOR (adding -lstdc++fs)..."
        sed -i 's/target_link_libraries(splitcode splitcode_core pthread)/target_link_libraries(splitcode splitcode_core pthread stdc++fs)/' src/CMakeLists.txt
    fi
fi
rm -rf build && mkdir -p build && cd build
cmake .. && make -j"$THREADS"
SPLITCODE_BIN="$WORKDIR/splitcode/build/src/splitcode"
echo "  splitcode binary: $SPLITCODE_BIN"

# ============================================================================
# 5. Set up analysis repo
# ============================================================================
echo "[5/8] Setting up analysis repo..."

# If this script is running from inside the cloned repo already, symlink it
# into WORKDIR so paths are consistent. Otherwise clone fresh.
if [ -f "$ANALYSIS_ROOT/scripts/run_all.sh" ]; then
    echo "  Running from cloned repo at $ANALYSIS_ROOT"
    if [ "$ANALYSIS_ROOT" != "$WORKDIR/seqproc-paper-analysis" ]; then
        ln -sfn "$ANALYSIS_ROOT" "$WORKDIR/seqproc-paper-analysis"
    fi
else
    if [ ! -d "$WORKDIR/seqproc-paper-analysis" ]; then
        git clone --branch phase3-orientation-benchmarks \
            git@github.com:COMBINE-lab/seqproc-paper-analysis.git \
            "$WORKDIR/seqproc-paper-analysis"
    fi
fi

cd "$WORKDIR/seqproc-paper-analysis"

# Python venv -- use the Python we identified in step 1
if [ -d "$MAMBA_ROOT_PREFIX/envs/bench" ] 2>/dev/null; then
    # micromamba env is already active, install directly
    pip install -q -r requirements.txt
else
    $PYTHON_BIN -m venv venv
    source venv/bin/activate
    pip install -q -r requirements.txt
fi

# ============================================================================
# 6. Download FULL SRA datasets
# ============================================================================
echo "[6/8] Downloading full SRA datasets..."
mkdir -p "$WORKDIR/seqproc-paper-analysis/data/10x_short"
cd "$WORKDIR/seqproc-paper-analysis/data"

# SPLiT-seq PE (SRR6750041) -- ~86.8M paired-end reads, ~20 GB
if [ ! -f SRR6750041_R1.fastq ]; then
    echo "  Downloading SRR6750041 (SPLiT-seq PE, ~20 GB)..."
    prefetch SRR6750041 && fasterq-dump --split-files SRR6750041 --threads "$THREADS"
    mv SRR6750041_1.fastq SRR6750041_R1.fastq
    mv SRR6750041_2.fastq SRR6750041_R2.fastq
    rm -rf SRR6750041/
    # Also create 10M subset used by some configs
    head -40000000 SRR6750041_R1.fastq > SRR6750041_10M_R1.fastq
    head -40000000 SRR6750041_R2.fastq > SRR6750041_10M_R2.fastq
fi

# LR-SPLiT-seq (SRR13948564) -- ~4.2M single-end long reads
if [ ! -f SRR13948564_full.fastq ]; then
    echo "  Downloading SRR13948564 (LR-SPLiT-seq, ~5 GB)..."
    prefetch SRR13948564 && fasterq-dump SRR13948564 --threads "$THREADS"
    mv SRR13948564.fastq SRR13948564_full.fastq
    rm -rf SRR13948564/
fi

# 10x Chromium v2 (SRR8315379) -- ~56.5M paired-end reads, ~10 GB
if [ ! -f 10x_short/SRR8315379_R1.fastq ]; then
    echo "  Downloading SRR8315379 (10x Chromium v2, ~10 GB)..."
    prefetch SRR8315379 && fasterq-dump --split-files SRR8315379 --threads "$THREADS"
    mv SRR8315379_1.fastq 10x_short/SRR8315379_R1.fastq
    mv SRR8315379_2.fastq 10x_short/SRR8315379_R2.fastq
    rm -rf SRR8315379/
fi

# sci-RNA-seq3 (SRR7827254) -- ~10.2M paired-end reads, ~3 GB
if [ ! -f SRR7827254_1.fastq ]; then
    echo "  Downloading SRR7827254 (sci-RNA-seq3, ~3 GB)..."
    prefetch SRR7827254 && fasterq-dump --split-files SRR7827254 --threads "$THREADS"
    # fasterq-dump names them _1.fastq / _2.fastq which matches data_config.py
    rm -rf SRR7827254/
fi

cd "$WORKDIR/seqproc-paper-analysis"

# Verify all data is present
echo "  Verifying data availability..."
python scripts/data_config.py --reads full

# ============================================================================
# 7. Run the full benchmark pipeline
# ============================================================================
echo "[7/8] Running full benchmark pipeline..."
export SEQPROC_BIN
export MATCHBOX_BIN
export SPLITCODE_BIN

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
