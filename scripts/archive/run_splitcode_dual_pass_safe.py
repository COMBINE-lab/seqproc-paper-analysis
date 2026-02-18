#!/usr/bin/env python3
"""
Dual-pass Splitcode benchmark script (Forward + Reverse Complement).
Uses file redirection to avoid memory issues and ensures robustness.
"""

import subprocess
import time
import os
import sys
from typing import Tuple, Set
from pathlib import Path

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
CONFIG_DIR = PROJECT_ROOT / "configs"

SPLITCODE_BIN = os.environ.get("SPLITCODE_BIN", str(PROJECT_ROOT.parent / "splitcode/build/src/splitcode"))
FULL_FQ = DATA_DIR / "SRR13948564_full.fastq"
FWD_CONFIG = CONFIG_DIR / "splitcode/splitseq_singleend_fwd.config"
NEW_CONFIG = CONFIG_DIR / "splitcode/splitseq_singleend_new.config"

def run_safe(cmd: str, log_file: str) -> Tuple[float, float, int]:
    """Run command with file redirection and return (runtime, peak_memory_mb, returncode)."""
    # Use /usr/bin/time -v to measure resources, redirecting ALL output to log_file
    full_cmd = f"/usr/bin/time -v {cmd} > {log_file} 2>&1"
    
    start = time.time()
    # We use check=False so we can handle errors manually
    result = subprocess.run(full_cmd, shell=True, cwd=PROJECT_ROOT)
    runtime = time.time() - start
    
    peak_mem_kb = 0
    if os.path.exists(log_file):
        with open(log_file, 'r', errors='replace') as f:
            for line in f:
                if 'Maximum resident set size' in line:
                    try:
                        peak_mem_kb = int(line.split(':')[1].strip())
                    except ValueError:
                        pass
    
    return runtime, peak_mem_kb / 1024, result.returncode

def get_valid_ids(fastq_path: str) -> Set[str]:
    """Extract valid read IDs from output FASTQ."""
    valid_ids = set()
    if not os.path.exists(fastq_path):
        return valid_ids
        
    try:
        with open(fastq_path, 'r') as f:
            while True:
                header = f.readline()
                if not header: break
                seq = f.readline()
                f.readline()
                f.readline()
                
                if header.startswith('@'):
                    # ID is first part of header after @
                    read_id = header.strip().split()[0][1:]
                    valid_ids.add(read_id)
    except Exception as e:
        print(f"Error reading {fastq_path}: {e}")
            
    return valid_ids

def main():
    print(f"Running Splitcode Dual-Pass Benchmark (Safe Mode)")
    print("=" * 60)
    
    # Check binaries and configs
    if not os.path.exists(SPLITCODE_BIN):
        print(f"Error: Splitcode binary not found at {SPLITCODE_BIN}")
        return
    if not FWD_CONFIG.exists():
        print(f"Error: Fwd config not found at {FWD_CONFIG}")
        return
    if not NEW_CONFIG.exists():
        print(f"Error: New config not found at {NEW_CONFIG}")
        return

    # Files
    fwd_out = "splitcode_fwd.fq"
    rc_out = "splitcode_rc.fq"
    mapping_fwd = "splitcode_fwd.map"
    mapping_rc = "splitcode_rc.map"
    log_fwd = "splitcode_fwd.log"
    log_rc = "splitcode_rc.log"

    # Cleanup
    for f in [fwd_out, rc_out, mapping_fwd, mapping_rc, log_fwd, log_rc]:
        if os.path.exists(f): os.remove(f)

    # 1. Forward Pass
    print(f"1. Running Forward Pass using {FWD_CONFIG}...")
    cmd_fwd = f"{SPLITCODE_BIN} -c {FWD_CONFIG} --assign -N 1 -t 8 -m {mapping_fwd} -o {fwd_out} {FULL_FQ}"
    rt_fwd, mem_fwd, rc_fwd = run_safe(cmd_fwd, log_fwd)
    
    print(f"   Time: {rt_fwd:.2f}s, Memory: {mem_fwd:.2f}MB, RC: {rc_fwd}")
    if rc_fwd != 0:
        print("   ERROR: Forward pass failed!")
        os.system(f"tail -n 20 {log_fwd}")
    elif os.path.exists(fwd_out):
        size_mb = os.path.getsize(fwd_out) / (1024*1024)
        print(f"   Output: {size_mb:.2f} MB")
        if size_mb < 10:
            print("   WARNING: Output too small!")
            os.system(f"head -n 20 {log_fwd}")

    # 2. Reverse Pass
    print(f"2. Running Reverse Pass using {NEW_CONFIG}...")
    cmd_rc = f"{SPLITCODE_BIN} -c {NEW_CONFIG} --assign -N 1 -t 8 -m {mapping_rc} -o {rc_out} {FULL_FQ}"
    rt_rc, mem_rc, rc_rc = run_safe(cmd_rc, log_rc)
    
    print(f"   Time: {rt_rc:.2f}s, Memory: {mem_rc:.2f}MB, RC: {rc_rc}")
    if rc_rc != 0:
        print("   ERROR: Reverse pass failed!")
        os.system(f"tail -n 20 {log_rc}")
    elif os.path.exists(rc_out):
        size_mb = os.path.getsize(rc_out) / (1024*1024)
        print(f"   Output: {size_mb:.2f} MB")
        if size_mb < 10:
            print("   WARNING: Output too small!")
            os.system(f"head -n 20 {log_rc}")
    
    # 3. Analyze Results
    print("3. Analyzing Results...")
    ids_fwd = get_valid_ids(fwd_out)
    ids_rc = get_valid_ids(rc_out)
    
    ids_union = ids_fwd.union(ids_rc)
    ids_intersect = ids_fwd.intersection(ids_rc)
    
    print("-" * 40)
    print(f"Forward Valid Reads: {len(ids_fwd):,}")
    print(f"Reverse Valid Reads: {len(ids_rc):,}")
    print(f"Union Valid Reads:   {len(ids_union):,}")
    print(f"Intersection:        {len(ids_intersect):,}")
    
    total_reads = 5_764_421
    recovery_pct = (len(ids_union) / total_reads) * 100
    
    print("-" * 40)
    print(f"Total Runtime: {rt_fwd + rt_rc:.2f}s")
    print(f"Peak Memory:   {max(mem_fwd, mem_rc):.2f}MB")
    print(f"Recovery:      {recovery_pct:.2f}%")

if __name__ == "__main__":
    main()
