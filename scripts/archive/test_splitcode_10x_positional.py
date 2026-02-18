#!/usr/bin/env python3
"""
Test Splitcode with positional extraction for 10x Chromium v2.
"""

import subprocess
import time
import os
import sys
from pathlib import Path

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "10x_short"
CONFIG_DIR = PROJECT_ROOT / "configs" / "splitcode"

SPLITCODE_BIN = os.environ.get("SPLITCODE_BIN", str(PROJECT_ROOT.parent / "splitcode/build/src/splitcode"))
R1_FQ = DATA_DIR / "SRR8315379_1M_R1.fastq"
CONFIG = CONFIG_DIR / "10x_v2_positional.config"

def run_with_memory(cmd: str) -> tuple:
    time_cmd = f"/usr/bin/time -v {cmd}"
    start = time.time()
    result = subprocess.run(time_cmd, shell=True, capture_output=True, text=True)
    runtime = time.time() - start
    
    peak_mem_kb = 0
    for line in result.stderr.split('\n'):
        if 'Maximum resident set size' in line:
            try:
                peak_mem_kb = int(line.split(':')[1].strip())
            except ValueError:
                pass
            break
            
    return runtime, peak_mem_kb / 1024, result.returncode, result.stderr, result.stdout

def main():
    print(f"Testing Splitcode 10x Positional Config")
    print(f"Config: {CONFIG}")
    print("=" * 60)
    
    out_fq = "splitcode_10x_pos.fq"
    if os.path.exists(out_fq): os.remove(out_fq)
    
    # Run splitcode
    # Note: 10x data is paired (R1 has tags, R2 has bio), but for benchmarking extraction we usually just check R1.
    # The user asked about "@extract 0:0<bc[16]><umi[10]>"
    cmd = f"{SPLITCODE_BIN} -c {CONFIG} -N 1 -t 4 -o {out_fq} {R1_FQ}"
    
    rt, mem, rc, err, out = run_with_memory(cmd)
    
    print(f"Runtime: {rt:.2f}s")
    print(f"Memory: {mem:.2f}MB")
    print(f"Return Code: {rc}")
    
    if rc != 0:
        print("Error output:")
        print(err)
        return
        
    if os.path.exists(out_fq):
        lines = 0
        with open(out_fq, 'rb') as f:
            for _ in f: lines += 1
        reads = lines // 4
        print(f"Reads extracted: {reads:,}")
        
        # Validation: Check first read structure
        with open(out_fq, 'r') as f:
            print("\nFirst read:")
            print(f.readline().strip())
            print(f.readline().strip())

if __name__ == "__main__":
    main()
