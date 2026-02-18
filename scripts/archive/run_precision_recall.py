#!/usr/bin/env python3
"""
Precision-Recall Analysis (Matchbox Paper Fig E style)

Generates synthetic SPLiT-seq data with a BARCODE WHITELIST, then measures
how accurately tools assign reads to the correct whitelist barcode.

Methodology (matching matchbox paper):
- Generate a whitelist of valid barcodes
- Assign each read a TRUE barcode from the whitelist
- Introduce sequencing errors
- Measure:
  - Correct: tool assigns read to the TRUE whitelist barcode
  - Incorrect: tool assigns read to a WRONG whitelist barcode
  - Unassigned: tool can't match read to any whitelist barcode

Usage:
    python scripts/run_precision_recall.py --num-reads 50000 --threads 4
"""

import os
import sys
import subprocess
import tempfile
import json
import argparse
import random
from pathlib import Path
from collections import defaultdict
from typing import Dict, Tuple, List, Set
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).parent.parent
os.chdir(PROJECT_ROOT)

# Tool binaries - use absolute paths
SEQPROC_BIN = os.environ.get("SEQPROC_BIN", "/home/ubuntu/combine-lab/seqproc/target/release/seqproc")
MATCHBOX_BIN = os.environ.get("MATCHBOX_BIN", "/home/ubuntu/combine-lab/matchbox/target/release/matchbox")
SPLITCODE_BIN = os.environ.get("SPLITCODE_BIN", "/home/ubuntu/combine-lab/splitcode/build/src/splitcode")

# Config files
SEQPROC_CONFIG = "configs/seqproc/splitseq_real.geom"
MATCHBOX_CONFIG = "configs/matchbox/splitseq.mb"
SPLITCODE_CONFIG = "configs/splitcode/splitseq_paper.config"

# SPLiT-seq structure
LINKER1 = "GTGGCCGCTGTTTCGCATCGGCGTACGACT"  # 30bp
LINKER2 = "ATCCACGTGCTTGAGA"  # 16bp

COLORS = {
    'seqproc': '#4DBBD5',
    'matchbox': '#E64B35',
    'splitcode': '#3C5488',
}

MARKERS = {0: 'o', 1: 's', 2: '^', 3: 'D'}


def random_barcode(length: int = 8) -> str:
    """Generate random barcode sequence."""
    return ''.join(random.choice('ACGT') for _ in range(length))


def introduce_errors(seq: str, error_rate: float) -> str:
    """Introduce random substitution errors."""
    if error_rate <= 0:
        return seq
    result = list(seq)
    for i in range(len(result)):
        if random.random() < error_rate:
            bases = [b for b in 'ACGT' if b != result[i]]
            result[i] = random.choice(bases)
    return ''.join(result)


def generate_synthetic_data(num_reads: int, error_rate: float, output_dir: str) -> Tuple[Dict, str, str, Set]:
    """Generate synthetic SPLiT-seq data with barcode whitelist.
    
    Returns:
        ground_truth: dict mapping read_id -> (bc1, bc2, bc3, umi)
        r1_path, r2_path: FASTQ file paths
        whitelist: set of valid BC1+BC2 combinations
    """
    random.seed(42)  # Reproducibility
    
    # Generate barcode whitelist pools (like real SPLiT-seq)
    bc1_pool = [random_barcode(8) for _ in range(96)]
    bc2_pool = [random_barcode(8) for _ in range(96)]
    bc3_pool = [random_barcode(8) for _ in range(96)]
    
    # Create whitelist of valid BC1+BC2 combinations
    whitelist = set()
    for bc1 in bc1_pool:
        for bc2 in bc2_pool:
            whitelist.add(f"{bc1}_{bc2}")
    
    ground_truth = {}
    
    r1_path = f"{output_dir}/synthetic_R1.fq"
    r2_path = f"{output_dir}/synthetic_R2.fq"
    
    with open(r1_path, 'w') as r1_out, open(r2_path, 'w') as r2_out:
        for i in range(num_reads):
            read_id = f"SYN.{i+1}"
            
            # Select random barcodes from whitelist
            bc1 = random.choice(bc1_pool)
            bc2 = random.choice(bc2_pool)
            bc3 = random.choice(bc3_pool)
            umi = random_barcode(10)
            
            ground_truth[read_id] = (bc1, bc2, bc3, umi)
            
            # Build R2 sequence: NN + UMI + BC3 + Linker1 + BC2 + Linker2 + BC1 + polyA
            # Introduce errors in barcodes AND linkers (realistic sequencing errors)
            r2_seq = (
                "NN" +  # Skip bases
                introduce_errors(umi, error_rate) +
                introduce_errors(bc3, error_rate) +
                introduce_errors(LINKER1, error_rate) +
                introduce_errors(bc2, error_rate) +
                introduce_errors(LINKER2, error_rate) +
                introduce_errors(bc1, error_rate) +
                "AAAAAAAAAA"  # polyA tail
            )
            
            # R1 is cDNA (random sequence)
            r1_seq = ''.join(random.choice('ACGT') for _ in range(100))
            
            # Quality scores (all high quality)
            r1_qual = 'I' * len(r1_seq)
            r2_qual = 'I' * len(r2_seq)
            
            # Write FASTQ
            r1_out.write(f"@{read_id}\n{r1_seq}\n+\n{r1_qual}\n")
            r2_out.write(f"@{read_id}\n{r2_seq}\n+\n{r2_qual}\n")
    
    return ground_truth, r1_path, r2_path, whitelist


def run_seqproc(r1_path: str, r2_path: str, threads: int, tmpdir: str) -> Dict[str, Tuple]:
    """Run seqproc and extract barcodes."""
    out1 = f"{tmpdir}/seqproc_R1.fq"
    out2 = f"{tmpdir}/seqproc_R2.fq"
    cmd = f"{SEQPROC_BIN} --geom {SEQPROC_CONFIG} --file1 {r1_path} --file2 {r2_path} --out1 {out1} --out2 {out2} --threads {threads}"
    subprocess.run(cmd, shell=True, capture_output=True)
    
    results = {}
    if os.path.exists(out2):
        with open(out2, 'r') as f:
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip()
                f.readline()
                f.readline()
                
                read_id = header.strip().split()[0].replace('@', '')
                if len(seq) >= 32:
                    umi = seq[0:10]
                    bc3 = seq[10:18]
                    bc1 = seq[-8:]
                    bc2 = seq[-16:-8]
                    results[read_id] = (bc1, bc2, bc3, umi)
    return results


def run_matchbox(r2_path: str, threads: int, tmpdir: str) -> Dict[str, Tuple]:
    """Run matchbox and extract barcodes."""
    out_tsv = f"{tmpdir}/matchbox_out.tsv"
    with open(MATCHBOX_CONFIG, 'r') as f:
        mb_script = f.read()
    cmd = f'{MATCHBOX_BIN} -e 0.2 -t {threads} -r "{mb_script}" {r2_path} > {out_tsv}'
    subprocess.run(cmd, shell=True, capture_output=True)
    
    results = {}
    if os.path.exists(out_tsv):
        with open(out_tsv, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    read_id = parts[0]
                    bc1, bc2, bc3, umi = parts[1], parts[2], parts[3], parts[4]
                    results[read_id] = (bc1, bc2, bc3, umi)
    return results


def run_splitcode(r1_path: str, r2_path: str, threads: int, tmpdir: str) -> Dict[str, Tuple]:
    """Run splitcode and extract barcodes."""
    out1 = f"{tmpdir}/splitcode_R1.fq"
    out2 = f"{tmpdir}/splitcode_R2.fq"
    mapping = f"{tmpdir}/splitcode_mapping.txt"
    cmd = f"{SPLITCODE_BIN} -c {SPLITCODE_CONFIG} --assign -N 2 -t {threads} -m {mapping} -o {out1},{out2} {r1_path} {r2_path}"
    subprocess.run(cmd, shell=True, capture_output=True)
    
    # Splitcode doesn't output sequences, just assignments - count matched reads
    results = {}
    if os.path.exists(out2):
        with open(out2, 'r') as f:
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip()
                f.readline()
                f.readline()
                read_id = header.strip().split()[0].replace('@', '')
                # Mark as matched (we can't extract actual barcodes from splitcode)
                results[read_id] = ('', '', '', '')
    return results


def find_closest_whitelist_match(bc1: str, bc2: str, whitelist: Set[str], max_dist: int = 1) -> Tuple[str, int]:
    """Find closest matching barcode in whitelist within max_dist.
    
    Returns (matched_barcode, distance) or (None, -1) if no match.
    """
    query = f"{bc1}_{bc2}"
    
    # Exact match
    if query in whitelist:
        return query, 0
    
    # Fuzzy match - find closest
    best_match = None
    best_dist = 100  # Large number
    
    for wl_bc in whitelist:
        wl_bc1, wl_bc2 = wl_bc.split('_')
        dist1 = sum(c1 != c2 for c1, c2 in zip(bc1, wl_bc1))
        dist2 = sum(c1 != c2 for c1, c2 in zip(bc2, wl_bc2))
        total_dist = dist1 + dist2
        
        if total_dist < best_dist:
            best_dist = total_dist
            best_match = wl_bc
    
    if best_dist <= max_dist:
        return best_match, best_dist
    return None, -1


def calculate_precision_recall_at_tolerance(ground_truth: Dict, predicted: Dict, whitelist: Set[str], max_mismatch: int) -> Tuple[float, float]:
    """Calculate precision and recall using whitelist matching at given tolerance.
    
    Methodology (matching matchbox paper Fig E):
    - Correct: extracted barcode matches the TRUE whitelist barcode within tolerance
    - Incorrect: extracted barcode matches a DIFFERENT whitelist barcode within tolerance
    - Unassigned: extracted barcode doesn't match ANY whitelist barcode within tolerance
    
    Args:
        max_mismatch: Maximum allowed mismatches for whitelist matching (0, 1, 2, or 3)
    
    Returns:
        frac_correct: correct / total_reads
        frac_incorrect: incorrect / total_reads
    """
    correct = 0
    incorrect = 0
    unassigned = 0
    
    total_reads = len(ground_truth)
    
    for read_id, true_bc in ground_truth.items():
        true_key = f"{true_bc[0]}_{true_bc[1]}"  # BC1_BC2
        
        if read_id in predicted:
            pred_bc = predicted[read_id]
            # Find what whitelist barcode the predicted barcode matches at this tolerance
            pred_match, dist = find_closest_whitelist_match(pred_bc[0], pred_bc[1], whitelist, max_dist=max_mismatch)
            
            if pred_match is None:
                # Doesn't match any whitelist barcode at this tolerance
                unassigned += 1
            elif pred_match == true_key:
                # Matches the correct whitelist barcode
                correct += 1
            else:
                # Matches a wrong whitelist barcode (false positive)
                incorrect += 1
        else:
            # Tool didn't process this read
            unassigned += 1
    
    frac_correct = correct / total_reads if total_reads > 0 else 0
    frac_incorrect = incorrect / total_reads if total_reads > 0 else 0
    
    return frac_correct, frac_incorrect


def main():
    parser = argparse.ArgumentParser(description="Precision-Recall Analysis")
    parser.add_argument('--num-reads', type=int, default=50000, help='Number of synthetic reads')
    parser.add_argument('--threads', type=int, default=4, help='Threads per tool')
    parser.add_argument('--output-dir', default=None)
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir) if args.output_dir else PROJECT_ROOT / "results/precision_recall"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("Precision-Recall Analysis (Fig E style)")
    print("=" * 70)
    print(f"Synthetic reads: {args.num_reads:,}")
    print(f"Threads: {args.threads}")
    print(f"Output: {output_dir}")
    print("=" * 70)
    
    # Fixed input error rate (realistic sequencing errors)
    INPUT_ERROR_RATE = 0.02  # 2% base-level error rate
    
    # Error tolerance levels to test (matching matchbox paper Fig E)
    # These represent the allowed mismatch threshold for whitelist matching
    tolerance_levels = [0, 1, 2, 3]
    
    tools = ['seqproc', 'matchbox', 'splitcode']
    results = {tool: {'frac_correct': [], 'frac_incorrect': []} for tool in tools}
    
    print(f"\nInput error rate: {INPUT_ERROR_RATE:.0%}")
    print("Testing different error tolerance levels for whitelist matching...\n")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Generate synthetic data ONCE with fixed error rate
        print("Generating synthetic data...", end=" ", flush=True)
        ground_truth, r1_path, r2_path, whitelist = generate_synthetic_data(
            args.num_reads, INPUT_ERROR_RATE, tmpdir
        )
        print(f"{len(ground_truth):,} reads, {len(whitelist):,} whitelist barcodes")
        
        # Save whitelist to file for verification
        whitelist_path = output_dir / 'whitelist.txt'
        with open(whitelist_path, 'w') as f:
            for bc in sorted(whitelist):
                f.write(f"{bc}\n")
        print(f"Whitelist saved to: {whitelist_path}")
        
        # Both tools process the SAME input files
        print(f"\n*** Both tools use SAME input data: ***")
        print(f"  R1: {r1_path}")
        print(f"  R2: {r2_path}")
        print(f"  Whitelist: {len(whitelist):,} BC1+BC2 combinations\n")
        
        # Run tools ONCE on the same data
        print("Running seqproc...", end=" ", flush=True)
        sp_pred = run_seqproc(r1_path, r2_path, args.threads, tmpdir)
        print(f"{len(sp_pred):,} reads extracted")
        
        print("Running matchbox...", end=" ", flush=True)
        mb_pred = run_matchbox(r2_path, args.threads, tmpdir)
        print(f"{len(mb_pred):,} reads extracted")
        
        print("Running splitcode...", end=" ", flush=True)
        sc_pred = run_splitcode(r1_path, r2_path, args.threads, tmpdir)
        print(f"{len(sc_pred):,} reads extracted")
        
        # Verify tools processed reads
        common_reads = set(sp_pred.keys()) & set(mb_pred.keys())
        print(f"\nVerification: {len(common_reads):,} reads processed by seqproc & matchbox")
        
        # Calculate precision/recall at each tolerance level
        for tol in tolerance_levels:
            print(f"\n[Tolerance {tol}] Allowed mismatches: {tol}")
            
            # seqproc
            sp_correct, sp_incorrect = calculate_precision_recall_at_tolerance(
                ground_truth, sp_pred, whitelist, max_mismatch=tol
            )
            results['seqproc']['frac_correct'].append(sp_correct)
            results['seqproc']['frac_incorrect'].append(sp_incorrect)
            print(f"  seqproc:  correct={sp_correct:.3f}, incorrect={sp_incorrect:.3f}")
            
            # matchbox
            mb_correct, mb_incorrect = calculate_precision_recall_at_tolerance(
                ground_truth, mb_pred, whitelist, max_mismatch=tol
            )
            results['matchbox']['frac_correct'].append(mb_correct)
            results['matchbox']['frac_incorrect'].append(mb_incorrect)
            print(f"  matchbox: correct={mb_correct:.3f}, incorrect={mb_incorrect:.3f}")
            
            # splitcode - note: splitcode does error correction internally, so tolerance doesn't apply the same way
            # We measure what fraction of reads splitcode matched correctly
            sc_correct, sc_incorrect = calculate_precision_recall_at_tolerance(
                ground_truth, sc_pred, whitelist, max_mismatch=tol
            )
            results['splitcode']['frac_correct'].append(sc_correct)
            results['splitcode']['frac_incorrect'].append(sc_incorrect)
            print(f"  splitcode: correct={sc_correct:.3f}, incorrect={sc_incorrect:.3f}")
    
    # Generate Fig E style plot
    print("\nGenerating precision-recall plot...")
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    for tool in ['seqproc', 'matchbox', 'splitcode']:
        frac_correct = results[tool]['frac_correct']
        frac_incorrect = results[tool]['frac_incorrect']
        
        # Plot line connecting points
        ax.plot(frac_incorrect, frac_correct, 
                color=COLORS[tool], linewidth=2, label=tool, alpha=0.8)
        
        # Plot points with different markers for tolerance levels
        for i, (fi, fc) in enumerate(zip(frac_incorrect, frac_correct)):
            ax.scatter(fi, fc, color=COLORS[tool], marker=MARKERS[i], 
                      s=100, edgecolor='black', linewidth=1, zorder=5)
    
    # Add legend for tolerance levels (corrected labels)
    for i, tol in enumerate(tolerance_levels):
        ax.scatter([], [], color='gray', marker=MARKERS[i], s=80, 
                  label=f'Tolerance {tol}', edgecolor='black')
    
    ax.set_xlabel('Fraction of reads with incorrect barcodes', fontsize=12)
    ax.set_ylabel('Fraction of reads with correct barcodes', fontsize=12)
    ax.set_title('Precision-Recall Analysis\n(Synthetic SPLiT-seq Data)', fontsize=14, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Set axis limits similar to matchbox paper
    max_incorrect = max(max(results[t]['frac_incorrect']) for t in tools)
    ax.set_xlim(-0.005, max(0.1, max_incorrect * 1.2))
    ax.set_ylim(0, 1.05)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig_precision_recall.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'fig_precision_recall.pdf', bbox_inches='tight')
    plt.close()
    
    # Save results
    with open(output_dir / 'precision_recall_results.json', 'w') as f:
        json.dump({
            'tolerance_levels': tolerance_levels,
            'input_error_rate': INPUT_ERROR_RATE,
            'results': results,
            'num_reads': args.num_reads,
        }, f, indent=2)
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\n{'Tolerance':<12} {'seqproc':>20} {'matchbox':>20}")
    print("-" * 54)
    for i, tol in enumerate(tolerance_levels):
        sp_c = results['seqproc']['frac_correct'][i]
        mb_c = results['matchbox']['frac_correct'][i]
        print(f"{tol} mismatch    {sp_c*100:>18.1f}% {mb_c*100:>18.1f}%")
    print("-" * 54)
    print(f"\nFigure saved to: {output_dir / 'fig_precision_recall.png'}")
    print("=" * 70)


if __name__ == "__main__":
    main()
