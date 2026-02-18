#!/bin/bash

# Script to generate whitelist files from barcode map TSVs.
# Run from the repository root: ./edit_distance_results/install_whitelists.sh

set -euo pipefail

echo "Generating whitelist files for seqproc..."

# Extract whitelists from map files into configs/seqproc/
echo "Creating splitseq_bc23_whitelist.txt from BC2 and BC3 maps..."
cut -f2 configs/seqproc/splitseq_bc2_seq2seq.tsv configs/seqproc/splitseq_bc3_seq2seq.tsv \
    > configs/seqproc/splitseq_bc23_whitelist.txt

echo "Creating splitseq_bc1_whitelist.txt from BC1 map..."
cut -f2 configs/seqproc/splitseq_bc1_seq2seq.tsv \
    > configs/seqproc/splitseq_bc1_whitelist.txt

echo "Whitelist files created at:"
echo "  configs/seqproc/splitseq_bc23_whitelist.txt"
echo "  configs/seqproc/splitseq_bc1_whitelist.txt"

# Verify
echo ""
echo "Verifying files:"
wc -l configs/seqproc/splitseq_bc23_whitelist.txt configs/seqproc/splitseq_bc1_whitelist.txt
