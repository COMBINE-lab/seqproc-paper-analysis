#!/bin/bash

# Script to install whitelist files at the expected locations
# This script should be run from the repository root

echo "Installing whitelist files for seqproc..."

# Create directory if it doesn't exist
sudo mkdir -p /home/ubuntu

# Extract whitelists from map files
echo "Creating whitelist_v1.txt from BC2 and BC3 maps..."
cut -f2 configs/seqproc/splitseq_bc2_seq2seq.tsv configs/seqproc/splitseq_bc3_seq2seq.tsv | cut -f2 | sudo tee /home/ubuntu/whitelist_v1.txt > /dev/null

echo "Creating bc1_whitelist.txt from BC1 map..."
cut -f2 configs/seqproc/splitseq_bc1_seq2seq.tsv | sudo tee /home/ubuntu/bc1_whitelist.txt > /dev/null

echo "Whitelist files installed at:"
echo "  /home/ubuntu/whitelist_v1.txt"
echo "  /home/ubuntu/bc1_whitelist.txt"

# Verify
echo -e "\nVerifying files:"
wc -l /home/ubuntu/whitelist_v1.txt /home/ubuntu/bc1_whitelist.txt
