#!/usr/bin/env python3
"""
SPLiT-seq Unified Benchmark + Accuracy Pipeline

Runs all three tools (seqproc, matchbox, splitcode) and generates:
1. Runtime benchmarks with box plots (configurable replicates)
2. Accuracy analysis with correlation plots
3. All publication-ready figures in one output directory

Usage:
    python scripts/run_splitseq_pipeline.py \
        --r1 /path/to/R1.fastq.gz \
        --r2 /path/to/R2.fastq.gz \
        --name SRR6750041_1M \
        --threads 4 \
        --replicates 5
"""

import os
import sys
import subprocess
import tempfile
import json
import argparse
import time
from pathlib import Path
from collections import defaultdict
from typing import Dict, Tuple, List, Optional
from dataclasses import dataclass, field
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).parent.parent
os.chdir(PROJECT_ROOT)

# Tool binaries (override with environment variables)
SEQPROC_BIN = os.environ.get("SEQPROC_BIN", "seqproc")
MATCHBOX_BIN = os.environ.get("MATCHBOX_BIN", "matchbox")
SPLITCODE_BIN = os.environ.get("SPLITCODE_BIN", "splitcode")

# Default config files (can be overridden via command line)
DEFAULT_SEQPROC_CONFIG = "configs/seqproc/splitseq_real.geom"
DEFAULT_MATCHBOX_CONFIG = "configs/matchbox/splitseq.mb"
DEFAULT_SPLITCODE_CONFIG = "configs/splitcode/splitseq_paper.config"

# Dataset presets for easy switching
DATASET_PRESETS = {
    'SRR6750041': {
        'description': 'Paired-end SPLiT-seq (77.6M reads)',
        'mode': 'paired',
        'seqproc_config': 'configs/seqproc/splitseq_real.geom',
        'matchbox_config': 'configs/matchbox/splitseq.mb',
        'splitcode_config': 'configs/splitcode/splitseq_paper.config',
    },
    'SRR13948564': {
        'description': 'Single-end SPLiT-seq (long-read)',
        'mode': 'single',
        'seqproc_config': 'configs/seqproc/splitseq_singleend.geom',
        'matchbox_config': 'configs/matchbox/splitseq_singleend.mb',
        'splitcode_config': 'configs/splitcode/splitseq_paper.config',
    },
    'ERR9958134': {
        'description': '10x Chromium GridION long-read',
        'mode': 'single',
        'seqproc_config': 'configs/seqproc/10x_longread.geom',
        'matchbox_config': 'configs/matchbox/10x_longread.mb',
        'splitcode_config': None,  # No splitcode config for 10x
    },
    'ERR9958135': {
        'description': '10x Chromium PromethION long-read',
        'mode': 'single',
        'seqproc_config': 'configs/seqproc/10x_longread.geom',
        'matchbox_config': 'configs/matchbox/10x_longread.mb',
        'splitcode_config': None,  # No splitcode config for 10x
    },
}

COLORS = {
    'seqproc': '#4DBBD5',
    'matchbox': '#E64B35',
    'splitcode': '#3C5488',
}


@dataclass
class PipelineConfig:
    name: str
    r1_path: str
    r2_path: str  # Same as r1_path for single-end
    threads: int
    replicates: int
    output_dir: Path
    seqproc_config: str = DEFAULT_SEQPROC_CONFIG
    matchbox_config: str = DEFAULT_MATCHBOX_CONFIG
    splitcode_config: str = DEFAULT_SPLITCODE_CONFIG
    mode: str = 'paired'  # 'paired' or 'single'
    total_reads: int = 0


@dataclass
class ToolResult:
    tool: str
    runtime: float
    reads_out: int
    returncode: int
    replicate: int


def count_fastq_reads(path: str) -> int:
    """Count reads in FASTQ file."""
    if not os.path.exists(path):
        return 0
    
    if path.endswith('.gz'):
        import gzip
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
    
    count = 0
    with opener(path, mode) as f:
        for _ in f:
            count += 1
    return count // 4


def run_command(cmd: str, cwd: Path = None) -> Tuple[float, int, str, str]:
    """Run command and return (runtime, returncode, stdout, stderr)."""
    start = time.time()
    result = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)
    runtime = time.time() - start
    return runtime, result.returncode, result.stdout, result.stderr


def run_seqproc(config: PipelineConfig, tmpdir: str) -> Tuple[float, int, str, str]:
    """Run seqproc."""
    out1 = f"{tmpdir}/seqproc_R1.fq"
    out2 = f"{tmpdir}/seqproc_R2.fq"
    
    if config.mode == 'single':
        # Native single-end mode: no --file2 or --out2
        cmd = f"{SEQPROC_BIN} --geom {config.seqproc_config} --file1 {config.r1_path} --out1 {out1} --threads {config.threads}"
        runtime, rc, stdout, stderr = run_command(cmd, PROJECT_ROOT)
        reads = count_fastq_reads(out1)
        return runtime, rc, reads, out1
    else:
        # Paired-end mode
        cmd = f"{SEQPROC_BIN} --geom {config.seqproc_config} --file1 {config.r1_path} --file2 {config.r2_path} --out1 {out1} --out2 {out2} --threads {config.threads}"
        runtime, rc, stdout, stderr = run_command(cmd, PROJECT_ROOT)
        reads = count_fastq_reads(out2)
        return runtime, rc, reads, out2


def run_matchbox(config: PipelineConfig, tmpdir: str) -> Tuple[float, int, str, str]:
    """Run matchbox using config file."""
    out_tsv = f"{tmpdir}/matchbox_out.tsv"
    # Read script from config file
    with open(config.matchbox_config, 'r') as f:
        mb_script = f.read()
    # For single-end, use r1_path (which equals r2_path)
    input_file = config.r1_path if config.mode == 'single' else config.r2_path
    cmd = f'{MATCHBOX_BIN} -e 0.2 -t {config.threads} -r "{mb_script}" {input_file} > {out_tsv}'
    runtime, rc, stdout, stderr = run_command(cmd, PROJECT_ROOT)
    
    # Count output lines
    reads = 0
    if os.path.exists(out_tsv):
        with open(out_tsv) as f:
            reads = sum(1 for _ in f)
    return runtime, rc, reads, out_tsv


def run_splitcode(config: PipelineConfig, tmpdir: str) -> Tuple[float, int, str, str]:
    """Run splitcode."""
    mapping = f"{tmpdir}/splitcode_mapping.txt"
    
    if config.mode == 'single':
        # Single-end mode
        out_fq = f"{tmpdir}/splitcode_out.fq"
        cmd = f"{SPLITCODE_BIN} -c {config.splitcode_config} --assign -N 1 -t {config.threads} -m {mapping} -o {out_fq} {config.r1_path}"
        runtime, rc, stdout, stderr = run_command(cmd, PROJECT_ROOT)
        reads = count_fastq_reads(out_fq)
        return runtime, rc, reads, out_fq
    else:
        # Paired-end mode
        out1 = f"{tmpdir}/splitcode_R1.fq"
        out2 = f"{tmpdir}/splitcode_R2.fq"
        cmd = f"{SPLITCODE_BIN} -c {config.splitcode_config} --assign -N 2 -t {config.threads} -m {mapping} -o {out1},{out2} {config.r1_path} {config.r2_path}"
        runtime, rc, stdout, stderr = run_command(cmd, PROJECT_ROOT)
        reads = count_fastq_reads(out2)
        return runtime, rc, reads, out2


def extract_seqproc_barcodes(fastq_path: str) -> Dict[str, Tuple[str, str, str, str]]:
    """Extract barcodes from seqproc output."""
    results = {}
    if not os.path.exists(fastq_path):
        return results
    with open(fastq_path, 'r') as f:
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


def extract_matchbox_barcodes(tsv_path: str) -> Dict[str, Tuple[str, str, str, str]]:
    """Extract barcodes from matchbox TSV output."""
    results = {}
    if not os.path.exists(tsv_path):
        return results
    with open(tsv_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                read_id = parts[0]
                bc1, bc2, bc3, umi = parts[1], parts[2], parts[3], parts[4]
                results[read_id] = (bc1, bc2, bc3, umi)
    return results


def extract_splitcode_barcodes(fastq_path: str) -> Dict[str, Tuple[str, str, str, str]]:
    """Extract barcodes from splitcode output.
    
    Splitcode assigns barcode IDs in the header, not sequences.
    We extract BC1+BC2 from the sequence for correlation analysis.
    """
    results = {}
    if not os.path.exists(fastq_path):
        return results
    with open(fastq_path, 'r') as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()
            f.readline()
            
            # Parse read ID from header
            read_id = header.strip().split()[0].replace('@', '')
            # Splitcode outputs the barcode-assigned sequence
            # Extract BC1+BC2 pattern similar to seqproc for correlation
            if len(seq) >= 16:
                bc1 = seq[-8:] if len(seq) >= 8 else seq
                bc2 = seq[-16:-8] if len(seq) >= 16 else ''
                bc3 = ''
                umi = ''
                results[read_id] = (bc1, bc2, bc3, umi)
    return results


def hamming_dist(s1: str, s2: str) -> int:
    """Calculate hamming distance."""
    if len(s1) != len(s2):
        return max(len(s1), len(s2))
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def calculate_agreement(sp_bc: Dict, mb_bc: Dict) -> Dict:
    """Calculate agreement between seqproc and matchbox."""
    common_reads = set(sp_bc.keys()) & set(mb_bc.keys())
    
    exact_match = 0
    bc_match_exact = 0
    bc_match_fuzzy = 0
    
    for read_id in common_reads:
        sp = sp_bc[read_id]
        mb = mb_bc[read_id]
        
        if sp == mb:
            exact_match += 1
        if sp[:3] == mb[:3]:
            bc_match_exact += 1
        
        bc1_ok = hamming_dist(sp[0], mb[0]) <= 1
        bc2_ok = hamming_dist(sp[1], mb[1]) <= 1
        bc3_ok = hamming_dist(sp[2], mb[2]) <= 1
        if bc1_ok and bc2_ok and bc3_ok:
            bc_match_fuzzy += 1
    
    return {
        'common_reads': len(common_reads),
        'exact_match': exact_match,
        'bc_match_exact': bc_match_exact,
        'bc_match_fuzzy': bc_match_fuzzy,
        'bc_match_fuzzy_rate': bc_match_fuzzy / len(common_reads) if common_reads else 0,
        'seqproc_only': len(set(sp_bc.keys()) - set(mb_bc.keys())),
        'matchbox_only': len(set(mb_bc.keys()) - set(sp_bc.keys())),
    }


def count_barcode_occurrences(barcodes: Dict) -> Dict[str, int]:
    """Count unique barcode combinations (BC1+BC2 only for cross-tool comparison)."""
    counts = defaultdict(int)
    for read_id, (bc1, bc2, bc3, umi) in barcodes.items():
        key = f"{bc1}_{bc2}"
        counts[key] += 1
    return dict(counts)


def generate_figures(results: List[ToolResult], accuracy: Dict, config: PipelineConfig):
    """Generate all publication-ready figures."""
    output_dir = config.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Determine which tools were run based on results
    tools_run = list(set(r.tool for r in results))
    tools = [t for t in ['seqproc', 'matchbox', 'splitcode'] if t in tools_run]
    
    # Aggregate results by tool
    tool_runtimes = {t: [] for t in tools}
    tool_reads = {t: 0 for t in tools}
    
    for r in results:
        tool_runtimes[r.tool].append(r.runtime)
        tool_reads[r.tool] = r.reads_out
    
    # === Figure 1: Runtime Box Plot ===
    fig, ax = plt.subplots(figsize=(8, 6))
    
    data = [tool_runtimes[t] for t in tools]
    bp = ax.boxplot(data, labels=[t.capitalize() for t in tools], patch_artist=True)
    
    for patch, tool in zip(bp['boxes'], tools):
        patch.set_facecolor(COLORS[tool])
        patch.set_alpha(0.7)
    
    ax.set_ylabel('Runtime (seconds)', fontsize=12)
    ax.set_title(f'SPLiT-seq Runtime Comparison\n{config.name} ({config.total_reads:,} reads, {config.replicates} reps)', 
                 fontsize=14, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    
    # Add mean values
    for i, tool in enumerate(tools):
        if tool_runtimes[tool]:
            mean_val = np.mean(tool_runtimes[tool])
            ax.text(i + 1, mean_val, f'{mean_val:.2f}s', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig_runtime_boxplot.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'fig_runtime_boxplot.pdf', bbox_inches='tight')
    plt.close()
    
    # === Figure 2: Speedup Bar Chart ===
    fig, ax = plt.subplots(figsize=(8, 6))
    
    mean_runtimes = {t: np.mean(tool_runtimes[t]) if tool_runtimes[t] else 0 for t in tools}
    baseline = mean_runtimes['matchbox'] if mean_runtimes['matchbox'] > 0 else 1
    speedups = [baseline / mean_runtimes[t] if mean_runtimes[t] > 0 else 0 for t in tools]
    
    bars = ax.bar([t.capitalize() for t in tools], speedups, color=[COLORS[t] for t in tools], edgecolor='black')
    
    for bar, speedup in zip(bars, speedups):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                f'{speedup:.2f}x', ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    ax.set_ylabel('Speedup (vs matchbox)', fontsize=12)
    ax.set_title(f'SPLiT-seq Speedup Comparison\n{config.name}', fontsize=14, fontweight='bold')
    ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig_speedup.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'fig_speedup.pdf', bbox_inches='tight')
    plt.close()
    
    # === Figure 3: Read Match Rate ===
    fig, ax = plt.subplots(figsize=(8, 6))
    
    match_rates = [tool_reads[t] / config.total_reads * 100 if config.total_reads > 0 else 0 for t in tools]
    bars = ax.bar([t.capitalize() for t in tools], match_rates, color=[COLORS[t] for t in tools], edgecolor='black')
    
    for bar, reads, rate in zip(bars, [tool_reads[t] for t in tools], match_rates):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{reads:,}\n({rate:.1f}%)', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    ax.set_ylabel('Read Match Rate (%)', fontsize=12)
    ax.set_title(f'SPLiT-seq Barcode Extraction\n{config.name} ({config.total_reads:,} total reads)', 
                 fontsize=14, fontweight='bold')
    ax.set_ylim(0, 100)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig_match_rate.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'fig_match_rate.pdf', bbox_inches='tight')
    plt.close()
    
    # === Figure 4: Barcode Correlation (seqproc vs matchbox) ===
    r_squared_mb = 0
    if 'sp_counts' in accuracy and 'mb_counts' in accuracy:
        sp_counts = accuracy['sp_counts']
        mb_counts = accuracy['mb_counts']
        common_bc = set(sp_counts.keys()) & set(mb_counts.keys())
        
        if len(common_bc) > 10:
            sp_vals = [sp_counts[bc] for bc in common_bc]
            mb_vals = [mb_counts[bc] for bc in common_bc]
            
            correlation = np.corrcoef(sp_vals, mb_vals)[0, 1]
            r_squared_mb = correlation ** 2
            
            fig, ax = plt.subplots(figsize=(8, 8))
            ax.scatter(mb_vals, sp_vals, alpha=0.4, s=15, c='#333333')
            
            max_val = max(max(sp_vals), max(mb_vals))
            ax.plot([1, max_val], [1, max_val], 'r--', alpha=0.7, linewidth=2, label='y=x')
            
            ax.set_xlabel('Reads per barcode (matchbox)', fontsize=12)
            ax.set_ylabel('Reads per barcode (seqproc)', fontsize=12)
            ax.set_title(f'Barcode Count Correlation (seqproc vs matchbox)\nR² = {r_squared_mb:.3f} ({len(common_bc):,} unique barcodes)', 
                         fontsize=14, fontweight='bold')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.legend(fontsize=11)
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(output_dir / 'fig_correlation_seqproc_matchbox.png', dpi=300, bbox_inches='tight')
            plt.savefig(output_dir / 'fig_correlation_seqproc_matchbox.pdf', bbox_inches='tight')
            plt.close()
    
    # === Figure 5: Barcode Correlation (seqproc vs splitcode) ===
    r_squared_sc = 0
    if 'sp_counts' in accuracy and 'sc_counts' in accuracy:
        sp_counts = accuracy['sp_counts']
        sc_counts = accuracy['sc_counts']
        common_bc = set(sp_counts.keys()) & set(sc_counts.keys())
        
        if len(common_bc) > 10:
            sp_vals = [sp_counts[bc] for bc in common_bc]
            sc_vals = [sc_counts[bc] for bc in common_bc]
            
            correlation = np.corrcoef(sp_vals, sc_vals)[0, 1]
            r_squared_sc = correlation ** 2
            
            fig, ax = plt.subplots(figsize=(8, 8))
            ax.scatter(sc_vals, sp_vals, alpha=0.4, s=15, c='#333333')
            
            max_val = max(max(sp_vals), max(sc_vals))
            ax.plot([1, max_val], [1, max_val], 'r--', alpha=0.7, linewidth=2, label='y=x')
            
            ax.set_xlabel('Reads per barcode (splitcode)', fontsize=12)
            ax.set_ylabel('Reads per barcode (seqproc)', fontsize=12)
            ax.set_title(f'Barcode Count Correlation (seqproc vs splitcode)\nR² = {r_squared_sc:.3f} ({len(common_bc):,} unique barcodes)', 
                         fontsize=14, fontweight='bold')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.legend(fontsize=11)
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(output_dir / 'fig_correlation_seqproc_splitcode.png', dpi=300, bbox_inches='tight')
            plt.savefig(output_dir / 'fig_correlation_seqproc_splitcode.pdf', bbox_inches='tight')
            plt.close()
    
    return r_squared_mb, r_squared_sc


def main():
    parser = argparse.ArgumentParser(
        description="SPLiT-seq Unified Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Dataset Presets:
  {', '.join(DATASET_PRESETS.keys())}

Examples:
  # Paired-end dataset (SRR6750041)
  python scripts/run_splitseq_pipeline.py --r1 R1.fq --r2 R2.fq --name SRR6750041 --preset SRR6750041
  
  # Single-end dataset (SRR13948564) 
  python scripts/run_splitseq_pipeline.py --r1 data.fq --name SRR13948564 --preset SRR13948564
"""
    )
    parser.add_argument('--r1', required=True, help='R1 FASTQ file (or single-end input)')
    parser.add_argument('--r2', default=None, help='R2 FASTQ file (optional for single-end)')
    parser.add_argument('--name', required=True, help='Dataset name')
    parser.add_argument('--preset', default=None, choices=list(DATASET_PRESETS.keys()),
                        help='Use preset configs for known datasets')
    parser.add_argument('--mode', default='paired', choices=['paired', 'single'],
                        help='Data mode: paired or single-end')
    parser.add_argument('--threads', type=int, default=4, help='Threads per tool')
    parser.add_argument('--replicates', type=int, default=5, help='Replicates for runtime benchmark')
    parser.add_argument('--seqproc-config', default=None, help='Override seqproc config')
    parser.add_argument('--matchbox-config', default=None, help='Override matchbox config')
    parser.add_argument('--splitcode-config', default=None, help='Override splitcode config')
    parser.add_argument('--output-dir', default=None)
    args = parser.parse_args()
    
    # Apply preset if specified
    if args.preset and args.preset in DATASET_PRESETS:
        preset = DATASET_PRESETS[args.preset]
        mode = preset['mode']
        seqproc_config = args.seqproc_config or preset['seqproc_config']
        matchbox_config = args.matchbox_config or preset['matchbox_config']
        splitcode_config = args.splitcode_config or preset['splitcode_config']
        print(f"Using preset: {args.preset} ({preset['description']})")
    else:
        mode = args.mode
        seqproc_config = args.seqproc_config or DEFAULT_SEQPROC_CONFIG
        matchbox_config = args.matchbox_config or DEFAULT_MATCHBOX_CONFIG
        splitcode_config = args.splitcode_config or DEFAULT_SPLITCODE_CONFIG
    
    # For single-end, r2 = r1
    r2_path = args.r2 if args.r2 else args.r1
    
    output_dir = Path(args.output_dir) if args.output_dir else PROJECT_ROOT / f"results/pipeline/{args.name}"
    
    # Count reads
    print("Counting input reads...")
    total_reads = count_fastq_reads(args.r1)
    
    config = PipelineConfig(
        name=args.name,
        r1_path=args.r1,
        r2_path=r2_path,
        threads=args.threads,
        replicates=args.replicates,
        output_dir=output_dir,
        seqproc_config=seqproc_config,
        matchbox_config=matchbox_config,
        splitcode_config=splitcode_config,
        mode=mode,
        total_reads=total_reads,
    )
    
    print("=" * 70)
    print("SPLiT-seq Unified Pipeline")
    print("=" * 70)
    print(f"Dataset:     {config.name}")
    print(f"Mode:        {config.mode}")
    print(f"R1:          {config.r1_path}")
    if config.mode == 'paired':
        print(f"R2:          {config.r2_path}")
    print(f"Total reads: {config.total_reads:,}")
    print(f"Threads:     {config.threads}")
    print(f"Replicates:  {config.replicates}")
    print(f"Configs:")
    print(f"  seqproc:   {config.seqproc_config}")
    print(f"  matchbox:  {config.matchbox_config}")
    print(f"  splitcode: {config.splitcode_config}")
    print(f"Output:      {config.output_dir}")
    print("=" * 70)
    
    all_results = []
    accuracy = {}
    
    # Run benchmarks with multiple replicates
    print(f"\n[1/2] Running benchmarks ({config.replicates} replicates)...")
    
    for rep in range(1, config.replicates + 1):
        print(f"\n  Replicate {rep}/{config.replicates}:")
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Seqproc
            print(f"    seqproc...", end=" ", flush=True)
            runtime, rc, reads, out_path = run_seqproc(config, tmpdir)
            print(f"{runtime:.2f}s, {reads:,} reads")
            all_results.append(ToolResult('seqproc', runtime, reads, rc, rep))
            
            # Extract barcodes on first rep for accuracy
            if rep == 1:
                accuracy['sp_barcodes'] = extract_seqproc_barcodes(out_path)
            
            # Matchbox
            print(f"    matchbox...", end=" ", flush=True)
            runtime, rc, reads, out_path = run_matchbox(config, tmpdir)
            print(f"{runtime:.2f}s, {reads:,} reads")
            all_results.append(ToolResult('matchbox', runtime, reads, rc, rep))
            
            if rep == 1:
                accuracy['mb_barcodes'] = extract_matchbox_barcodes(out_path)
            
            # Splitcode (optional - skip if no config)
            if config.splitcode_config:
                print(f"    splitcode...", end=" ", flush=True)
                runtime, rc, reads, out_path = run_splitcode(config, tmpdir)
                print(f"{runtime:.2f}s, {reads:,} reads")
                all_results.append(ToolResult('splitcode', runtime, reads, rc, rep))
                
                if rep == 1:
                    accuracy['sc_barcodes'] = extract_splitcode_barcodes(out_path)
    
    # Calculate accuracy metrics
    print("\n[2/2] Calculating accuracy...")
    
    sp_bc = accuracy.get('sp_barcodes', {})
    mb_bc = accuracy.get('mb_barcodes', {})
    sc_bc = accuracy.get('sc_barcodes', {})
    
    agreement = calculate_agreement(sp_bc, mb_bc)
    accuracy['agreement'] = agreement
    accuracy['sp_counts'] = count_barcode_occurrences(sp_bc)
    accuracy['mb_counts'] = count_barcode_occurrences(mb_bc)
    accuracy['sc_counts'] = count_barcode_occurrences(sc_bc)
    
    print(f"  Common reads (seqproc/matchbox): {agreement['common_reads']:,}")
    print(f"  Fuzzy BC match: {agreement['bc_match_fuzzy']:,} ({agreement['bc_match_fuzzy_rate']:.1%})")
    
    # Generate figures
    print("\nGenerating figures...")
    r_squared_mb, r_squared_sc = generate_figures(all_results, accuracy, config)
    
    # Save summary
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Aggregate stats
    tools = ['seqproc', 'matchbox', 'splitcode']
    summary = {
        'config': {
            'name': config.name,
            'total_reads': config.total_reads,
            'threads': config.threads,
            'replicates': config.replicates,
        },
        'runtime': {},
        'accuracy': {
            'seqproc_reads': len(sp_bc),
            'matchbox_reads': len(mb_bc),
            'splitcode_reads': ([r.reads_out for r in all_results if r.tool == 'splitcode'] or [0])[0],
            'common_reads': agreement['common_reads'],
            'bc_match_fuzzy_rate': agreement['bc_match_fuzzy_rate'],
            'correlation_r_squared_mb': r_squared_mb,
            'correlation_r_squared_sc': r_squared_sc,
        }
    }
    
    for tool in tools:
        runtimes = [r.runtime for r in all_results if r.tool == tool]
        reads = [r.reads_out for r in all_results if r.tool == tool]
        summary['runtime'][tool] = {
            'mean': np.mean(runtimes) if runtimes else 0,
            'std': np.std(runtimes) if runtimes else 0,
            'reads': reads[0] if reads else 0,
        }
    
    with open(output_dir / 'pipeline_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Save commands
    with open(output_dir / 'commands.txt', 'w') as f:
        f.write(f"# SPLiT-seq Pipeline Commands\n# Mode: {config.mode}\n\n")
        if config.mode == 'single':
            f.write(f"# seqproc (single-end)\n{SEQPROC_BIN} --geom {config.seqproc_config} --file1 {config.r1_path} --out1 OUT_R1.fq --out2 OUT_R2.fq --threads {config.threads}\n\n")
            f.write(f"# matchbox\n{MATCHBOX_BIN} -e 0.2 -t {config.threads} -r '...see config...' {config.r1_path}\n\n")
            f.write(f"# splitcode (single-end)\n{SPLITCODE_BIN} -c {config.splitcode_config} --assign -N 1 -t {config.threads} -m mapping.txt -o OUT.fq {config.r1_path}\n")
        else:
            f.write(f"# seqproc\n{SEQPROC_BIN} --geom {config.seqproc_config} --file1 {config.r1_path} --file2 {config.r2_path} --out1 OUT_R1.fq --out2 OUT_R2.fq --threads {config.threads}\n\n")
            f.write(f"# matchbox\n{MATCHBOX_BIN} -e 0.2 -t {config.threads} -r '...see config...' {config.r2_path}\n\n")
            f.write(f"# splitcode\n{SPLITCODE_BIN} -c {config.splitcode_config} --assign -N 2 -t {config.threads} -m mapping.txt -o OUT_R1.fq,OUT_R2.fq {config.r1_path} {config.r2_path}\n")
    
    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    print(f"\n{'Tool':<12} {'Runtime':>10} {'Std':>8} {'Reads':>12} {'Match%':>8} {'Speedup':>8}")
    print("-" * 60)
    
    baseline = summary['runtime']['matchbox']['mean']
    for tool in tools:
        r = summary['runtime'][tool]
        match_pct = r['reads'] / config.total_reads * 100 if config.total_reads > 0 else 0
        speedup = baseline / r['mean'] if r['mean'] > 0 else 0
        print(f"{tool:<12} {r['mean']:>9.2f}s {r['std']:>7.2f}s {r['reads']:>12,} {match_pct:>7.1f}% {speedup:>7.2f}x")
    
    print("-" * 60)
    print(f"Barcode correlation (seqproc vs matchbox): R² = {r_squared_mb:.4f}")
    if 'splitcode' in tools:
        print(f"Barcode correlation (seqproc vs splitcode): R² = {r_squared_sc:.4f}")
    print(f"\nFigures saved to: {output_dir}")
    print("=" * 70)


if __name__ == "__main__":
    main()
