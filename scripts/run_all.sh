#!/bin/bash
#
# Master benchmark script for seqproc paper analysis.
#
# Runs the complete pipeline: performance benchmarks, concordance analysis,
# discordant read validation, and figure generation.
#
# Usage:
#   ./scripts/run_all.sh              # 1M-read subsets (default)
#   ./scripts/run_all.sh --reads full # Full SRA datasets
#   ./scripts/run_all.sh --reads 1m --threads 8
#
# Prerequisites:
#   - Tool binaries: set SEQPROC_BIN, MATCHBOX_BIN, SPLITCODE_BIN env vars
#     (or ensure they are at the default paths in ../combine-lab/seqproc/ etc.)
#   - Python environment: pip install -r requirements.txt
#   - FASTQ data files in data/ (run: python scripts/data_config.py --reads <level>
#     to check availability)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Defaults
READS="1m"
THREADS=4
REPLICATES=3

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --reads)
            READS="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --replicates)
            REPLICATES="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [--reads 1m|full] [--threads N] [--replicates N]"
            echo ""
            echo "Options:"
            echo "  --reads      Dataset size: '1m' (default) or 'full'"
            echo "  --threads    Number of threads (default: 4)"
            echo "  --replicates Number of benchmark replicates (default: 3)"
            echo ""
            echo "Steps executed:"
            echo "  1. Check data availability"
            echo "  2. Run performance benchmarks (Table 2)"
            echo "  3. Run concordance analysis (Figures 3-5)"
            echo "  4. Run discordant read validation (Supp. Figure S1)"
            echo "  5. Generate all publication figures"
            exit 0
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Run with --help for usage."
            exit 1
            ;;
    esac
done

# Determine Python interpreter (use whatever python3 is on PATH;
# setup_and_run.sh ensures the micromamba bench env is on PATH)
PYTHON="python3"

echo "========================================================================"
echo "SEQPROC PAPER ANALYSIS -- FULL PIPELINE"
echo "========================================================================"
echo "  Reads level:  $READS"
echo "  Threads:      $THREADS"
echo "  Replicates:   $REPLICATES"
echo "  Python:       $PYTHON"
echo "  Project root: $PROJECT_ROOT"
echo ""

# Step 0: Check data availability
echo "--- Step 0: Checking data availability ---"
"$PYTHON" "$SCRIPT_DIR/data_config.py" --reads "$READS"
echo ""

# Step 1: Performance benchmarks (Table 2)
echo "========================================================================"
echo "--- Step 1: Performance Benchmarks (Table 2) ---"
echo "========================================================================"
"$PYTHON" "$SCRIPT_DIR/run_paper_benchmarks.py" \
    --threads "$THREADS" \
    --replicates "$REPLICATES" \
    --reads "$READS" \
    --datasets splitseq_pe_raw splitseq_se_raw 10x_short sciseq

# Step 2: Concordance analysis (Figures 3-5 data)
echo ""
echo "========================================================================"
echo "--- Step 2: Concordance Analysis (Figures 3-5) ---"
echo "========================================================================"
"$PYTHON" "$SCRIPT_DIR/concordance_analysis.py" \
    --threads "$THREADS" \
    --reads "$READS"

# Step 3: Discordant read structural validation (Supp. Figure S1)
echo ""
echo "========================================================================"
echo "--- Step 3: Discordant Read Validation (Supp. Fig S1) ---"
echo "========================================================================"
"$PYTHON" "$SCRIPT_DIR/discordant_analysis.py" --reads "$READS"

# Step 4: Generate all publication figures
echo ""
echo "========================================================================"
echo "--- Step 4: Generate Publication Figures ---"
echo "========================================================================"
"$PYTHON" "$SCRIPT_DIR/generate_figures.py"

echo ""
echo "========================================================================"
echo "PIPELINE COMPLETE"
echo "========================================================================"
echo "Results:  $PROJECT_ROOT/results/paper_figures/"
echo "Figures:  fig_concordance_heatmaps.pdf"
echo "          fig_recovery_comparison.pdf"
echo "          fig_hamming_vs_edit.pdf"
echo "          fig_discordant_summary.pdf"
echo "JSON:     benchmark_results.json"
echo ""
echo "To update the paper figures:"
echo "  cp results/paper_figures/*.pdf ../paper/Figures/"
