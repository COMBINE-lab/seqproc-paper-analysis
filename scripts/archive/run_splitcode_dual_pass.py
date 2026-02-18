#!/usr/bin/env python3
"""
Dual-pass Splitcode benchmark script (Forward + Reverse Complement).
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

def get_valid_ids(fastq_path: str) -> Set[str]:
    """Extract valid read IDs from output FASTQ."""
    valid_ids = set()
    if not os.path.exists(fastq_path):
        return valid_ids
        
    with open(fastq_path, 'r') as f:
        while True:
            header = f.readline()
            if not header: break
            seq = f.readline()
            f.readline()
            f.readline()
            
            # ID is first part of header after @
            # e.g. @SRR13948564.2 -> SRR13948564.2
            read_id = header.strip().split()[0][1:]
            valid_ids.add(read_id)
            
    return valid_ids

def main():
    print(f"Running Splitcode Dual-Pass Benchmark")
    print("=" * 60)
    
    # Files
    fwd_out = "splitcode_fwd.fq"
    rc_out = "splitcode_rc.fq" # The "new" config is for RC (or whatever the second pass is)
    # Actually, the user said "those splitcode numbers are only for reverse" for the previous run (using new config).
    # So "new config" = Reverse/RC pass
    # "fwd config" = Forward pass (splitseq_singleend_fwd.config)
    
    mapping_fwd = "splitcode_fwd.map"
    mapping_rc = "splitcode_rc.map"

    # Cleanup
    for f in [fwd_out, rc_out, mapping_fwd, mapping_rc]:
        p = PROJECT_ROOT / f
        if p.exists(): p.unlink()

    # Use absolute paths for outputs
    fwd_out_abs = PROJECT_ROOT / fwd_out
    rc_out_abs = PROJECT_ROOT / rc_out
    mapping_fwd_abs = PROJECT_ROOT / mapping_fwd
    mapping_rc_abs = PROJECT_ROOT / mapping_rc

    # Verify configuration files exist
    if not FWD_CONFIG.exists():
        print(f"Error: Forward configuration file '{FWD_CONFIG}' does not exist.")
        sys.exit(1)
    if not NEW_CONFIG.exists():
        print(f"Error: Reverse configuration file '{NEW_CONFIG}' does not exist.")
        sys.exit(1)

    # 1. Forward Pass
    print(f"1. Running Forward Pass using {FWD_CONFIG}...")
    cmd_fwd = f"{SPLITCODE_BIN} -c {FWD_CONFIG} --assign -N 1 -t 8 -m {mapping_fwd} -o {fwd_out} {FULL_FQ}"
    rt_fwd, mem_fwd, rc_fwd, err_fwd = run_with_memory(cmd_fwd, PROJECT_ROOT)
    print(f"   Time: {rt_fwd:.2f}s, Memory: {mem_fwd:.2f}MB, RC: {rc_fwd}")
    if rc_fwd != 0 or os.path.getsize(fwd_out) < 1000000:
        print("   WARNING: Forward pass output suspicious!")
        print(f"   Stderr: {err_fwd[:1000]}...")

    # 2. Reverse Pass (using 'new' config provided by user)
    print(f"2. Running Reverse Pass using {NEW_CONFIG}...")
    cmd_rc = f"{SPLITCODE_BIN} -c {NEW_CONFIG} --assign -N 1 -t 8 -m {mapping_rc} -o {rc_out} {FULL_FQ}"
    rt_rc, mem_rc, rc_rc, err_rc = run_with_memory(cmd_rc, PROJECT_ROOT)
    print(f"   Time: {rt_rc:.2f}s, Memory: {mem_rc:.2f}MB, RC: {rc_rc}")
    if rc_rc != 0 or os.path.getsize(rc_out) < 1000000:
        print("   WARNING: Reverse pass output suspicious!")
        print(f"   Stderr: {err_rc[:1000]}...")

    # Check output size
    if fwd_out_abs.exists():
        print(f"   Output size: {fwd_out_abs.stat().st_size / (1024*1024):.2f} MB")
    if rc_out_abs.exists():
        print(f"   Output size: {rc_out_abs.stat().st_size / (1024*1024):.2f} MB")

    # 3. Analyze Results
    print("3. Analyzing Results...")
    ids_fwd = get_valid_ids(str(fwd_out_abs))
    ids_rc = get_valid_ids(str(rc_out_abs))
    
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