#!/usr/bin/env python3
"""
Generate barcode correlation figure for FULL dataset (seqproc vs matchbox).
"""

import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from pathlib import Path

# Paths to full dataset outputs
SEQPROC_FQ = "/home/ubuntu/seqproc_full_fwd.fq"
MATCHBOX_TSV = "/home/ubuntu/matchbox_full.tsv"
OUTPUT_DIR = Path("results/paper_figures_full")

def parse_seqproc(fq_path):
    print(f"Parsing Seqproc output: {fq_path} ...")
    barcodes = {}
    count = 0
    with open(fq_path, 'r') as f:
        while True:
            header = f.readline()
            if not header: break
            seq = f.readline().strip()
            f.readline() # +
            f.readline() # qual
            
            # Seqproc Geom: 1{<umi><bc3><bc2><bc1>}
            # Lengths: 10 + 8 + 8 + 8 = 34
            # We compare BC1 and BC2 as per original script
            if len(seq) >= 26:
                # BC1 is last 8
                bc1 = seq[-8:]
                # BC2 is 8 before that
                bc2 = seq[-16:-8]
                
                # Read ID: @SRR... 123/1 -> SRR...
                # Original script: read_id = header.strip().split()[0].replace('@', '')
                read_id = header.strip().split()[0].replace('@', '')
                barcodes[read_id] = f"{bc1}_{bc2}"
                count += 1
                if count % 500000 == 0:
                    print(f"  Parsed {count} reads...")
    print(f"Total Seqproc reads parsed: {len(barcodes):,}")
    return barcodes

def parse_matchbox(tsv_path):
    print(f"Parsing Matchbox output: {tsv_path} ...")
    barcodes = {}
    count = 0
    with open(tsv_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                # Format: read_id \t bc1 \t bc2 \t bc3 \t umi
                read_id = parts[0]
                bc1 = parts[1]
                bc2 = parts[2]
                
                barcodes[read_id] = f"{bc1}_{bc2}"
                count += 1
                if count % 500000 == 0:
                    print(f"  Parsed {count} reads...")
    print(f"Total Matchbox reads parsed: {len(barcodes):,}")
    return barcodes

def main():
    sp_barcodes = parse_seqproc(SEQPROC_FQ)
    mb_barcodes = parse_matchbox(MATCHBOX_TSV)
    
    print("Counting barcodes...")
    sp_counts = defaultdict(int)
    for bc in sp_barcodes.values():
        sp_counts[bc] += 1
        
    mb_counts = defaultdict(int)
    for bc in mb_barcodes.values():
        mb_counts[bc] += 1
        
    # Find common barcodes
    common_bc = set(sp_counts.keys()) & set(mb_counts.keys())
    print(f"Common barcodes: {len(common_bc):,}")
    
    if len(common_bc) > 10:
        sp_vals = [sp_counts[bc] for bc in common_bc]
        mb_vals = [mb_counts[bc] for bc in common_bc]
        
        correlation = np.corrcoef(sp_vals, mb_vals)[0, 1]
        r_squared = correlation ** 2
        
        print(f"R-squared: {r_squared}")
        
        fig, ax = plt.subplots(figsize=(8, 8))
        # Plot with transparency to show density
        ax.scatter(mb_vals, sp_vals, alpha=0.4, s=15, c='#333333')
        
        max_val = max(max(sp_vals), max(mb_vals))
        ax.plot([1, max_val], [1, max_val], 'r--', alpha=0.7, linewidth=2, label='y=x')
        
        ax.set_xlabel('Reads per barcode (matchbox)', fontsize=12)
        ax.set_ylabel('Reads per barcode (seqproc)', fontsize=12)
        ax.set_title(f'Barcode Count Correlation (seqproc vs matchbox)\nR² = {r_squared:.3f} ({len(common_bc):,} unique barcodes)', 
                     fontsize=14, fontweight='bold')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        plt.tight_layout()
        out_png = OUTPUT_DIR / 'fig7_barcode_correlation.png'
        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {out_png}")
    else:
        print("Not enough common barcodes!")

if __name__ == "__main__":
    main()
