#!/bin/bash
set -euo pipefail

###############################################################################
# seqproc Full Paper Benchmark -- Clean Machine Setup + Run
#
# Assumes: Ubuntu Linux, internet access, sudo privileges.
# Produces: All paper figures and benchmark data for FULL SRA datasets.
#
# Usage: scp this to the target machine, then:
#   chmod +x setup_and_run.sh && ./setup_and_run.sh
###############################################################################

WORKDIR="$HOME/seqproc-bench"
THREADS=$(nproc)
REPLICATES=3

echo "================================================================"
echo "seqproc Full Paper Benchmark Setup"
echo "  WORKDIR:    $WORKDIR"
echo "  THREADS:    $THREADS"
echo "  REPLICATES: $REPLICATES"
echo "================================================================"

# ============================================================================
# 1. System dependencies
# ============================================================================
echo "[1/8] Installing system dependencies..."
sudo apt-get update -qq
sudo apt-get install -y -qq \
    build-essential cmake curl git python3 python3-pip python3-venv \
    time pkg-config libssl-dev

# Rust toolchain
if ! command -v cargo &>/dev/null; then
    echo "  Installing Rust..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    source "$HOME/.cargo/env"
fi
# Ensure cargo is on PATH for this session
export PATH="$HOME/.cargo/bin:$PATH"

# SRA toolkit (fasterq-dump)
if ! command -v fasterq-dump &>/dev/null; then
    echo "  Installing SRA toolkit..."
    curl -L -o /tmp/sratoolkit.tar.gz \
        https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
    tar -xzf /tmp/sratoolkit.tar.gz -C /tmp
    SRA_DIR=$(ls -d /tmp/sratoolkit.* | head -1)
    sudo cp "$SRA_DIR/bin/fasterq-dump" /usr/local/bin/
    sudo cp "$SRA_DIR/bin/prefetch" /usr/local/bin/
    rm -rf /tmp/sratoolkit*
fi

mkdir -p "$WORKDIR"

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
mkdir -p build && cd build
cmake .. && make -j"$THREADS"
SPLITCODE_BIN="$WORKDIR/splitcode/build/src/splitcode"
echo "  splitcode binary: $SPLITCODE_BIN"

# ============================================================================
# 5. Clone paper analysis repo
# ============================================================================
echo "[5/8] Cloning paper analysis repo..."
if [ ! -d "$WORKDIR/seqproc-paper-analysis" ]; then
    git clone --branch phase3-orientation-benchmarks \
        git@github.com:COMBINE-lab/seqproc-paper-analysis.git \
        "$WORKDIR/seqproc-paper-analysis"
fi

cd "$WORKDIR/seqproc-paper-analysis"

# Python venv
python3 -m venv venv
source venv/bin/activate
pip install -q -r requirements.txt

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
