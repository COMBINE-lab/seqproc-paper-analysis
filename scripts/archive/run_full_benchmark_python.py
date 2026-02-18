#!/usr/bin/env python3
"""
Benchmark script for Full Dataset (SRR13948564) - Fwd and RC passes.
Uses the same timing/memory measurement logic as run_paper_benchmarks.py.
"""

import subprocess
import time
import os
from typing import Tuple
from pathlib import Path

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
CONFIG_DIR = PROJECT_ROOT / "configs"

SEQPROC_BIN = os.environ.get("SEQPROC_BIN", str(PROJECT_ROOT.parent / "seqproc/target/release/seqproc"))
MATCHBOX_BIN = os.environ.get("MATCHBOX_BIN", str(PROJECT_ROOT.parent / "matchbox/target/release/matchbox"))

FULL_FQ = DATA_DIR / "SRR13948564_full.fastq"

def run_with_memory(cmd: str, cwd=None) -> Tuple[float, float, int, str]:
    """Run command and return (runtime, peak_memory_mb, returncode, stderr)."""
    time_cmd = f"/usr/bin/time -v {cmd}"
    start = time.time()
    result = subprocess.run(time_cmd, shell=True, capture_output=True, text=True, cwd=cwd)
    runtime = time.time() - start
    
    peak_mem_kb = 0
    for line in result.stderr.split('\n'):
        if 'Maximum resident set size' in line:
            try:
                peak_mem_kb = int(line.split(':')[1].strip())
            except ValueError:
                pass
            break
    
    return runtime, peak_mem_kb / 1024, result.returncode, result.stderr

def main():
    print(f"Benchmarking Full Dataset: {FULL_FQ}")
    print("=" * 60)

    n_replicates = 3
    results = {
        'seqproc': {'time': [], 'mem': []},
        'matchbox': {'time': [], 'mem': []}
    }

    for i in range(n_replicates):
        print(f"\nReplicate {i+1}/{n_replicates}")
        print("-" * 20)

        # 1. Seqproc Forward
        print("Running Seqproc Forward...")
        cmd_sp_fwd = f"{SEQPROC_BIN} --geom {CONFIG_DIR}/seqproc/splitseq_singleend_primer.geom --file1 {FULL_FQ} --out1 /dev/null --threads 8"
        rt_sp_fwd, mem_sp_fwd, rc_sp_fwd, _ = run_with_memory(cmd_sp_fwd, PROJECT_ROOT)
        
        # 2. Seqproc RC
        print("Running Seqproc RC...")
        cmd_sp_rc = f"{SEQPROC_BIN} --geom {CONFIG_DIR}/seqproc/splitseq_singleend_rc.geom --file1 {FULL_FQ} --out1 /dev/null --threads 8"
        rt_sp_rc, mem_sp_rc, rc_sp_rc, _ = run_with_memory(cmd_sp_rc, PROJECT_ROOT)
        
        sp_total_time = rt_sp_fwd + rt_sp_rc
        sp_max_mem = max(mem_sp_fwd, mem_sp_rc)
        results['seqproc']['time'].append(sp_total_time)
        results['seqproc']['mem'].append(sp_max_mem)
        print(f"  Seqproc Dual: {sp_total_time:.2f}s, {sp_max_mem:.2f}MB")

        # 3. Matchbox Forward
        print("Running Matchbox Forward...")
        cmd_mb_fwd = f"{MATCHBOX_BIN} -e 0.2 -t 8 -s {CONFIG_DIR}/matchbox/splitseq_singleend.mb {FULL_FQ} > /dev/null"
        rt_mb_fwd, mem_mb_fwd, rc_mb_fwd, _ = run_with_memory(cmd_mb_fwd, PROJECT_ROOT)
        
        # 4. Matchbox RC
        print("Running Matchbox RC...")
        cmd_mb_rc = f"{MATCHBOX_BIN} -e 0.2 -t 8 -s {CONFIG_DIR}/matchbox/splitseq_singleend_rc.mb {FULL_FQ} > /dev/null"
        rt_mb_rc, mem_mb_rc, rc_mb_rc, _ = run_with_memory(cmd_mb_rc, PROJECT_ROOT)
        
        mb_total_time = rt_mb_fwd + rt_mb_rc
        mb_max_mem = max(mem_mb_fwd, mem_mb_rc)
        results['matchbox']['time'].append(mb_total_time)
        results['matchbox']['mem'].append(mb_max_mem)
        print(f"  Matchbox Dual: {mb_total_time:.2f}s, {mb_max_mem:.2f}MB")

    print("\n" + "=" * 60)
    print(f"Summary ({n_replicates} Replicates):")
    
    sp_avg_time = sum(results['seqproc']['time']) / n_replicates
    sp_avg_mem = sum(results['seqproc']['mem']) / n_replicates
    print(f"Seqproc:  Avg Time = {sp_avg_time:.2f}s, Avg Memory = {sp_avg_mem:.2f}MB")
    
    mb_avg_time = sum(results['matchbox']['time']) / n_replicates
    mb_avg_mem = sum(results['matchbox']['mem']) / n_replicates
    print(f"Matchbox: Avg Time = {mb_avg_time:.2f}s, Avg Memory = {mb_avg_mem:.2f}MB")

if __name__ == "__main__":
    main()
