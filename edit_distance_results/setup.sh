#!/bin/bash

# Setup script for edit distance benchmark
# This script helps set up the environment for running the benchmark

echo "=== Edit Distance Benchmark Setup ==="
echo

# Check if we're in the right directory
if [ ! -f "run_paper_benchmarks.py" ]; then
    echo "Error: Please run this script from the edit_distance_results directory"
    exit 1
fi

# Set up environment variables
echo "Setting up environment variables..."

# Default paths - adjust if your tools are elsewhere
export SEQPROC_BIN="${SEQPROC_BIN:-$(cd ../.. && pwd)/combine-lab/seqproc/target/release/seqproc}"
export MATCHBOX_BIN="${MATCHBOX_BIN:-$(cd ../.. && pwd)/matchbox/target/release/matchbox}"
export SPLITCODE_BIN="${SPLITCODE_BIN:-$(cd ../.. && pwd)/splitcode/build/src/splitcode}"

echo "SEQPROC_BIN: $SEQPROC_BIN"
echo "MATCHBOX_BIN: $MATCHBOX_BIN"
echo "SPLITCODE_BIN: $SPLITCODE_BIN"
echo

# Check if binaries exist
echo "Checking binaries..."
for bin in "$SEQPROC_BIN" "$MATCHBOX_BIN" "$SPLITCODE_BIN"; do
    if [ -f "$bin" ]; then
        echo "✓ Found: $bin"
    else
        echo "✗ Missing: $bin"
        echo "  Please install the tool or set the environment variable"
    fi
done
echo

# Check data directory
DATA_DIR="$(cd ../.. && pwd)/data"
echo "Checking data directory: $DATA_DIR"

if [ -d "$DATA_DIR" ]; then
    echo "✓ Data directory exists"
    
    # Check for required files
    echo "Checking required data files..."
    required_files=(
        "data/SRR6750041_1M_R1.fastq"
        "data/SRR6750041_1M_R2.fastq"
        "data/SRR13948564_1M.fastq"
        "data/10x_short/SRR8315379_1M_R1.fastq"
        "data/10x_short/SRR8315379_1M_R2.fastq"
        "data/SRR7827254_1M_1.fastq"
        "data/SRR7827254_1M_2.fastq"
        "data/10x/ERR9958134_1M.fastq"
        "data/10x/ERR9958135_1M.fastq"
        "data/3M-february-2018.txt.gz"
    )
    
    for file in "${required_files[@]}"; do
        if [ -f "$DATA_DIR/../$file" ]; then
            echo "✓ Found: $file"
        else
            echo "✗ Missing: $file"
        fi
    done
else
    echo "✗ Data directory not found"
fi
echo

# Check Python dependencies
echo "Checking Python dependencies..."
python3 -c "import numpy, matplotlib; print('✓ Python dependencies satisfied')" 2>/dev/null || {
    echo "✗ Python dependencies missing"
    echo "  Install with: sudo apt install python3-numpy python3-matplotlib"
}
echo

# Check whitelist files
echo "Checking whitelist files..."
if [ -f "configs/whitelist_v1.txt" ] && [ -f "configs/bc1_whitelist.txt" ]; then
    echo "✓ Whitelist files found"
else
    echo "✗ Whitelist files missing"
    echo "  Creating from map files..."
    
    # Create whitelists if they don't exist
    if [ -f "../configs/seqproc/splitseq_bc1_seq2seq.tsv" ]; then
        cut -f2 ../configs/seqproc/splitseq_bc1_seq2seq.tsv > configs/bc1_whitelist.txt
        cut -f2 ../configs/seqproc/splitseq_bc2_seq2seq.tsv ../configs/seqproc/splitseq_bc3_seq2seq.tsv | cut -f2 > configs/whitelist_v1.txt
        echo "✓ Whitelist files created"
    fi
fi
echo

echo "=== Setup Complete ==="
echo
echo "To run the benchmark:"
echo "  cd $(pwd)"
echo "  python3 run_paper_benchmarks.py --threads 4 --replicates 1"
echo
echo "For full results (3 replicates):"
echo "  python3 run_paper_benchmarks.py --threads 4 --replicates 3"
