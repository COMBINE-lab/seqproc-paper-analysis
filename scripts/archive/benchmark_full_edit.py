#!/usr/bin/env python3
import subprocess
import time
import os
import argparse
from pathlib import Path
from typing import Dict, List, Tuple

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent
SEQPROC_BIN = "/home/ubuntu/combine-lab/seqproc/target/release/seqproc"
DATA_DIR = PROJECT_ROOT / "data"
CONFIG_DIR = PROJECT_ROOT / "configs/seqproc"

# Dataset Definitions
DATASETS = {
    'splitseq_se': {
        'name': 'SPLiT-seq SE Long (SRR13948564)',
        'file': DATA_DIR / "SRR13948564_1M.fastq",
        'mode': 'single',
        'configs': {
            'Hamming': CONFIG_DIR / "splitseq_singleend_primer.geom",
            'Edit': CONFIG_DIR / "splitseq_singleend_primer_edit.geom"
        }
    },
    '10x_gridion': {
        'name': '10x GridION (ERR9958134)',
        'file': DATA_DIR / "10x/ERR9958134_1M.fastq",
        'mode': 'dual_pass', # Fwd + Rev
        'configs': {
            'Hamming': {
                'fwd': CONFIG_DIR / "10x_longread_fwd.geom",
                'rev': CONFIG_DIR / "10x_longread_rev.geom"
            },
            'Edit': {
                'fwd': CONFIG_DIR / "10x_longread_fwd_edit.geom",
                'rev': CONFIG_DIR / "10x_longread_rev_edit.geom"
            }
        }
    },
    '10x_promethion': {
        'name': '10x PromethION (ERR9958135)',
        'file': DATA_DIR / "10x/ERR9958135_1M.fastq",
        'mode': 'dual_pass',
        'configs': {
            'Hamming': {
                'fwd': CONFIG_DIR / "10x_longread_fwd.geom",
                'rev': CONFIG_DIR / "10x_longread_rev.geom"
            },
            'Edit': {
                'fwd': CONFIG_DIR / "10x_longread_fwd_edit.geom",
                'rev': CONFIG_DIR / "10x_longread_rev_edit.geom"
            }
        }
    },
    'splitseq_pe_replace': {
        'name': 'SPLiT-seq PE Replace (SRR6750041)',
        'r1': DATA_DIR / "SRR6750041_1M_R1.fastq",
        'r2': DATA_DIR / "SRR6750041_1M_R2.fastq",
        'mode': 'paired',
        'configs': {
            'Hamming': CONFIG_DIR / "splitseq_replacement.geom",
            'Edit': CONFIG_DIR / "splitseq_replacement_edit.geom"
        },
        'maps': [
            CONFIG_DIR / "splitseq_bc3_seq2seq.tsv",
            CONFIG_DIR / "splitseq_bc2_seq2seq.tsv",
            CONFIG_DIR / "splitseq_bc1_seq2seq.tsv"
        ]
    },
    'sciseq3': {
        'name': 'Sci-Seq 3 (SRR7827254)',
        'r1': DATA_DIR / "SRR7827254_1M_1.fastq",
        'r2': DATA_DIR / "SRR7827254_1M_2.fastq",
        'mode': 'paired',
        'configs': {
            'Hamming': CONFIG_DIR / "sciseq3.geom",
            'Edit': CONFIG_DIR / "sciseq3_edit.geom"
        }
    }
}

def run_cmd(cmd):
    start = time.time()
    try:
        subprocess.run(cmd, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e.stderr.decode()}")
        return 0, 0
    return time.time() - start, 0

def count_reads(filepath):
    if not os.path.exists(filepath):
        return 0
    try:
        proc = subprocess.run(["wc", "-l", filepath], capture_output=True, text=True)
        lines = int(proc.stdout.split()[0])
        return lines // 4
    except:
        return 0

def run_benchmark(dataset_key, threads=4):
    ds = DATASETS[dataset_key]
    print(f"\nBenchmarking {ds['name']}...")
    
    results = {}
    
    for method in ['Hamming', 'Edit']:
        print(f"  Running {method}...", end=" ", flush=True)
        
        config = ds['configs'][method]
        runtime = 0
        reads = 0
        
        if ds['mode'] == 'single':
            output_file = f"bench_{dataset_key}_{method.lower()}.fq"
            cmd = [
                SEQPROC_BIN,
                "-g", str(config),
                "-1", str(ds['file']),
                "-o", output_file,
                "-t", str(threads)
            ]
            rt, _ = run_cmd(cmd)
            runtime = rt
            reads = count_reads(output_file)
            if os.path.exists(output_file): os.remove(output_file)
            
        elif ds['mode'] == 'paired':
            out1 = f"bench_{dataset_key}_{method.lower()}_R1.fq"
            out2 = f"bench_{dataset_key}_{method.lower()}_R2.fq"
            cmd = [
                SEQPROC_BIN,
                "-g", str(config),
                "-1", str(ds['r1']),
                "-2", str(ds['r2']),
                "-o", out1,
                "-w", out2, # Using -w for out2 based on help message
                "-t", str(threads)
            ]
            
            # Add map files if present
            if 'maps' in ds:
                for m in ds['maps']:
                    cmd.extend(["-a", str(m)])
            
            rt, _ = run_cmd(cmd)
            runtime = rt
            # Usually we count valid output reads (R2 usually contains barcodes for SPLiT-seq/SciSeq?)
            # Just sum them or pick one? Let's pick R1 for simplicity as both should have same count if paired
            reads = count_reads(out1)
            
            if os.path.exists(out1): os.remove(out1)
            if os.path.exists(out2): os.remove(out2)
            
        elif ds['mode'] == 'dual_pass':
            # Run Fwd
            out_fwd_1 = f"bench_{dataset_key}_{method.lower()}_fwd_R1.fq"
            cmd_fwd = [
                SEQPROC_BIN,
                "-g", str(config['fwd']),
                "-1", str(ds['file']),
                "-o", out_fwd_1,
                "-t", str(threads)
            ]
            rt1, _ = run_cmd(cmd_fwd)
            
            # Run Rev
            out_rev_1 = f"bench_{dataset_key}_{method.lower()}_rev_R1.fq"
            cmd_rev = [
                SEQPROC_BIN,
                "-g", str(config['rev']),
                "-1", str(ds['file']),
                "-o", out_rev_1,
                "-t", str(threads)
            ]
            rt2, _ = run_cmd(cmd_rev)
            
            runtime = rt1 + rt2
            # Combine reads (count R1s)
            reads = count_reads(out_fwd_1) + count_reads(out_rev_1)
            
            # Cleanup parts
            for f in [out_fwd_1, out_rev_1]:
                if os.path.exists(f): os.remove(f)
            
        print(f"Done. ({runtime:.2f}s, {reads:,} reads)")
        results[method] = {'time': runtime, 'reads': reads}

    # Report
    h = results['Hamming']
    e = results['Edit']
    
    delta_reads = e['reads'] - h['reads']
    pct_reads = (delta_reads / h['reads'] * 100) if h['reads'] > 0 else 0
    delta_time = e['time'] - h['time']
    speedup = h['time'] / e['time'] if e['time'] > 0 else 0
    
    print("-" * 60)
    print(f"{'Metric':<15} | {'Hamming':<10} | {'Edit':<10} | {'Delta':<10}")
    print("-" * 60)
    print(f"{'Time (s)':<15} | {h['time']:<10.2f} | {e['time']:<10.2f} | {speedup:.2f}x speedup")
    print(f"{'Reads':<15} | {h['reads']:<10,} | {e['reads']:<10,} | {pct_reads:+.1f}%")
    print("-" * 60)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    
    for key in DATASETS:
        run_benchmark(key, args.threads)

if __name__ == "__main__":
    main()
