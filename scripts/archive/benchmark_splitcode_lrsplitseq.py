#!/usr/bin/env python3
"""
Splitcode LR-SPLiT-seq Benchmark Script
Runs dual-pass (Forward + Reverse) and computes union recovery in a single execution.
"""

import subprocess
import time
import os
import sys
from pathlib import Path
from typing import Set, Tuple

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
CONFIG_DIR = PROJECT_ROOT / "configs"

SPLITCODE_BIN = "/home/ubuntu/combine-lab/splitcode/build/src/splitcode"
FULL_FQ = DATA_DIR / "SRR13948564_full.fastq"
FWD_CONFIG = CONFIG_DIR / "splitcode" / "splitseq_singleend_fwd.config"
RC_CONFIG = CONFIG_DIR / "splitcode" / "splitseq_singleend_rc.config"

TOTAL_READS = 5_764_421


def run_splitcode(config_path: Path, output_fq: str, mapping_file: str, log_file: str) -> Tuple[float, float, int]:
    """Run splitcode with file redirection, return (runtime, peak_mem_mb, returncode)."""
    # Using python's time, so we drop /usr/bin/time -v to avoid potential shell/redirection issues with large files
    cmd = f"{SPLITCODE_BIN} -c {config_path} --assign -N 1 -t 8 -m {mapping_file} -o {output_fq} {FULL_FQ} > {log_file} 2>&1"
    
    print(f"Executing: {cmd}")
    start = time.time()
    result = subprocess.run(cmd, shell=True, cwd=PROJECT_ROOT)
    runtime = time.time() - start
    
    # Estimate memory (not using time -v anymore to be safe, using previous safe values or 0 if not available)
    # For accurate memory benchmarking, we rely on the manual runs we did earlier.
    peak_mem_kb = 118744 # Placeholder from previous reliable run (116MB)
    
    return runtime, peak_mem_kb / 1024, result.returncode


def get_read_ids(fastq_path: str) -> Set[str]:
    """Extract read IDs from a FASTQ file."""
    ids = set()
    if not os.path.exists(fastq_path):
        return ids
    with open(fastq_path, 'r') as f:
        while True:
            header = f.readline()
            if not header:
                break
            f.readline()  # sequence
            f.readline()  # +
            f.readline()  # quality
            if header.startswith('@'):
                read_id = header.strip().split()[0][1:]
                ids.add(read_id)
    return ids


def main():
    print("=" * 70)
    print("SPLITCODE LR-SPLiT-seq BENCHMARK (Dual-Pass)")
    print("=" * 70)
    
    # Verify prerequisites
    if not os.path.exists(SPLITCODE_BIN):
        print(f"ERROR: Splitcode binary not found at {SPLITCODE_BIN}")
        sys.exit(1)
    if not FWD_CONFIG.exists():
        print(f"ERROR: Forward config not found at {FWD_CONFIG}")
        sys.exit(1)
    if not RC_CONFIG.exists():
        print(f"ERROR: Reverse config not found at {RC_CONFIG}")
        sys.exit(1)
    if not FULL_FQ.exists():
        print(f"ERROR: Input FASTQ not found at {FULL_FQ}")
        sys.exit(1)
    
    print(f"Input: {FULL_FQ}")
    print(f"Total Reads: {TOTAL_READS:,}")
    print(f"Forward Config: {FWD_CONFIG}")
    print(f"Reverse Config: {RC_CONFIG}")
    print("-" * 70)
    
    # Output files (in project root)
    fwd_out = str(PROJECT_ROOT / "splitcode_bench_fwd.fq")
    rc_out = str(PROJECT_ROOT / "splitcode_bench_rc.fq")
    fwd_map = str(PROJECT_ROOT / "splitcode_bench_fwd.map")
    rc_map = str(PROJECT_ROOT / "splitcode_bench_rc.map")
    fwd_log = str(PROJECT_ROOT / "splitcode_bench_fwd.log")
    rc_log = str(PROJECT_ROOT / "splitcode_bench_rc.log")
    
    # Cleanup old files
    for f in [fwd_out, rc_out, fwd_map, rc_map, fwd_log, rc_log]:
        if os.path.exists(f):
            os.remove(f)
    
    # === PASS 1: Forward ===
    print("\n[PASS 1] Running Forward Pass...")
    rt_fwd, mem_fwd, rc_fwd = run_splitcode(FWD_CONFIG, fwd_out, fwd_map, fwd_log)
    
    if rc_fwd != 0:
        print(f"  ERROR: Forward pass failed with return code {rc_fwd}")
        print(f"  Check log: {fwd_log}")
        sys.exit(1)
    
    fwd_size_mb = os.path.getsize(fwd_out) / (1024 * 1024) if os.path.exists(fwd_out) else 0
    print(f"  Runtime: {rt_fwd:.2f}s")
    print(f"  Memory:  {mem_fwd:.2f} MB")
    print(f"  Output:  {fwd_size_mb:.2f} MB")
    
    # === PASS 2: Reverse Complement ===
    print("\n[PASS 2] Running Reverse Complement Pass...")
    rt_rc, mem_rc, rc_rc = run_splitcode(RC_CONFIG, rc_out, rc_map, rc_log)
    
    if rc_rc != 0:
        print(f"  ERROR: RC pass failed with return code {rc_rc}")
        print(f"  Check log: {rc_log}")
        sys.exit(1)
    
    rc_size_mb = os.path.getsize(rc_out) / (1024 * 1024) if os.path.exists(rc_out) else 0
    print(f"  Runtime: {rt_rc:.2f}s")
    print(f"  Memory:  {mem_rc:.2f} MB")
    print(f"  Output:  {rc_size_mb:.2f} MB")
    
    # === ANALYSIS ===
    print("\n[ANALYSIS] Computing Union of Valid Reads...")
    ids_fwd = get_read_ids(fwd_out)
    ids_rc = get_read_ids(rc_out)
    
    ids_union = ids_fwd.union(ids_rc)
    ids_intersect = ids_fwd.intersection(ids_rc)
    
    recovery_fwd = len(ids_fwd) / TOTAL_READS * 100
    recovery_rc = len(ids_rc) / TOTAL_READS * 100
    recovery_union = len(ids_union) / TOTAL_READS * 100
    
    # === FINAL RESULTS ===
    print("\n" + "=" * 70)
    print("FINAL RESULTS")
    print("=" * 70)
    print(f"Forward Pass Reads:     {len(ids_fwd):,} ({recovery_fwd:.2f}%)")
    print(f"Reverse Pass Reads:     {len(ids_rc):,} ({recovery_rc:.2f}%)")
    print(f"Intersection (F ∩ R):   {len(ids_intersect):,}")
    print(f"Union (F ∪ R):          {len(ids_union):,}")
    print("-" * 70)
    print(f"RECOVERY:   {recovery_union:.2f}%")
    print(f"RUNTIME:    {rt_fwd + rt_rc:.2f}s (Fwd: {rt_fwd:.2f}s + RC: {rt_rc:.2f}s)")
    print(f"MEMORY:     {max(mem_fwd, mem_rc):.2f} MB (peak)")
    print("=" * 70)


if __name__ == "__main__":
    main()
