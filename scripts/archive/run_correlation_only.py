#!/usr/bin/env python3
"""
Run only the barcode correlation analysis - to resume after crash.
Uses subsetting to avoid memory issues.
"""

import os
import sys
import subprocess
import tempfile
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
os.chdir(PROJECT_ROOT)

SEQPROC_BIN = os.environ.get('SEQPROC_BIN', '../seqproc/target/release/seqproc')
MATCHBOX_BIN = os.environ.get('MATCHBOX_BIN', '../matchbox/target/release/matchbox')

# Use 1M subset for correlation to avoid memory crash
DATASET = {
    'r1': 'data/SRR6750041_1M_R1.fastq',
    'r2': 'data/SRR6750041_1M_R2.fastq',
    'seqproc_geom': 'configs/seqproc/splitseq_real.geom',
    'matchbox_config': 'configs/matchbox/splitseq.mb',
}

OUTPUT_DIR = Path('results/paper_figures_full')


def run_barcode_correlation_analysis(threads: int = 4):
    """Generate barcode count correlation plots for SPLiT-seq data."""
    
    print("="*70)
    print("Running Barcode Correlation Analysis (1M subset for memory safety)")
    print("="*70)
    
    if not os.path.exists(DATASET['r1']) or not os.path.exists(DATASET['r2']):
        print(f"Data not found: {DATASET['r1']}")
        return
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Run seqproc
        print("Running seqproc for barcode extraction...", end=" ", flush=True)
        sp_out1 = f"{tmpdir}/seqproc_R1.fq"
        sp_out2 = f"{tmpdir}/seqproc_R2.fq"
        cmd = f"{SEQPROC_BIN} --geom {DATASET['seqproc_geom']} --file1 {DATASET['r1']} --file2 {DATASET['r2']} --out1 {sp_out1} --out2 {sp_out2} --threads {threads}"
        subprocess.run(cmd, shell=True, capture_output=True, cwd=PROJECT_ROOT)
        
        sp_barcodes = {}
        if os.path.exists(sp_out2):
            with open(sp_out2, 'r') as f:
                while True:
                    header = f.readline()
                    if not header:
                        break
                    seq = f.readline().strip()
                    f.readline()
                    f.readline()
                    read_id = header.strip().split()[0].replace('@', '')
                    if len(seq) >= 32:
                        bc1 = seq[-8:]
                        bc2 = seq[-16:-8]
                        sp_barcodes[read_id] = f"{bc1}_{bc2}"
        print(f"{len(sp_barcodes):,} reads")
        
        # Run matchbox
        print("Running matchbox for barcode extraction...", end=" ", flush=True)
        mb_out = f"{tmpdir}/matchbox_out.tsv"
        with open(DATASET['matchbox_config'], 'r') as f:
            mb_script = f.read()
        cmd = f'{MATCHBOX_BIN} -e 0.2 -t {threads} -r "{mb_script}" {DATASET["r2"]} > {mb_out}'
        subprocess.run(cmd, shell=True, capture_output=True, cwd=PROJECT_ROOT)
        
        mb_barcodes = {}
        if os.path.exists(mb_out):
            with open(mb_out, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        read_id = parts[0]
                        bc1, bc2 = parts[1], parts[2]
                        mb_barcodes[read_id] = f"{bc1}_{bc2}"
        print(f"{len(mb_barcodes):,} reads")
        
        # Count barcode occurrences
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
            
            fig, ax = plt.subplots(figsize=(8, 8))
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
            plt.savefig(OUTPUT_DIR / 'fig7_barcode_correlation.png', dpi=300, bbox_inches='tight')
            plt.savefig(OUTPUT_DIR / 'fig7_barcode_correlation.pdf', bbox_inches='tight')
            plt.close()
            
            print(f"\nCorrelation figure saved: R² = {r_squared:.3f}")
            print(f"Output: {OUTPUT_DIR / 'fig7_barcode_correlation.png'}")
        else:
            print("Not enough common barcodes for correlation plot")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    
    run_barcode_correlation_analysis(args.threads)
    print("\nDone!")
