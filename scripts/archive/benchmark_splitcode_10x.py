#!/usr/bin/env python3
"""
Benchmark Splitcode on 10x v2 Short Reads (1M subset)
"""
import subprocess
import time
import os
import sys
from pathlib import Path

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent
SPLITCODE_BIN = "/home/ubuntu/combine-lab/splitcode/build/src/splitcode"
CONFIG_FILE = PROJECT_ROOT / "configs/splitcode/10x_v2_user.config"
INPUT_FQ = PROJECT_ROOT / "data/10x_short/SRR8315379_1M_R1.fastq"
OUTPUT_FQ = PROJECT_ROOT / "splitcode_10x_final.out"
LOG_FILE = PROJECT_ROOT / "splitcode_10x_benchmark.log"

def main():
    print(f"Benchmarking Splitcode on {INPUT_FQ}...")
    
    # Clean previous run
    if OUTPUT_FQ.exists(): os.remove(OUTPUT_FQ)
    if LOG_FILE.exists(): os.remove(LOG_FILE)

    cmd = f"/usr/bin/time -v {SPLITCODE_BIN} --x-only -c {CONFIG_FILE} --pipe {INPUT_FQ} > {OUTPUT_FQ} 2> {LOG_FILE}"
    
    start = time.time()
    result = subprocess.run(cmd, shell=True)
    runtime = time.time() - start
    
    # Parse memory
    peak_mem_kb = 0
    if LOG_FILE.exists():
        with open(LOG_FILE) as f:
            for line in f:
                if "Maximum resident set size" in line:
                    try:
                        peak_mem_kb = int(line.split(':')[1].strip())
                    except: pass
    
    # Check lines
    line_count = 0
    if OUTPUT_FQ.exists():
        # FASTQ is 4 lines per read. wc -l gives total lines. 
        # But wait, splitcode output here might be FASTQ. 
        # The user said "seems to output the matched bc and umi stiched together (i.e. the input read)"
        # Let's count lines.
        import subprocess as sp
        out = sp.check_output(f"wc -l {OUTPUT_FQ}", shell=True).decode().strip().split()[0]
        line_count = int(out)
        
    reads = line_count / 4
    
    print(f"Runtime: {runtime:.2f}s")
    print(f"Memory: {peak_mem_kb/1024:.2f} MB")
    print(f"Reads Output: {reads}")
    print(f"Recovery: {(reads/1000000)*100:.2f}%")

if __name__ == "__main__":
    main()
