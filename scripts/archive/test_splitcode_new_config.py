#!/usr/bin/env python3
"""
Test script for Splitcode with new configuration on LR-SPLiT-seq full dataset.
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
NEW_CONFIG = CONFIG_DIR / "splitcode/splitseq_singleend_new.config"

# Maps for validation
BC1_MAP = CONFIG_DIR / "seqproc/splitseq_bc1_seq2seq.tsv"
BC2_MAP = CONFIG_DIR / "seqproc/splitseq_bc2_seq2seq.tsv"
BC3_MAP = CONFIG_DIR / "seqproc/splitseq_bc3_seq2seq.tsv"

class SplitSeqSingleEndValidityAnalyzer:
    """Analyzes SPLiT-seq Single-End reads for validity (d<=1)."""
    
    LINKER1 = "GTGGCCGATGTTTCGCATCGGCGTACGACT"  # 30bp
    LINKER2 = "ATCCACGTGCTTGAGACTGTGG" 
    
    def __init__(self, bc1_map, bc2_map, bc3_map):
        self.bc1_wl = self._load_whitelist(bc1_map)
        self.bc2_wl = self._load_whitelist(bc2_map)
        self.bc3_wl = self._load_whitelist(bc3_map)
        
    def _load_whitelist(self, path):
        wl = set()
        with open(path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    wl.add(parts[1])
        return wl
        
    def _hamming(self, s1, s2):
        if len(s1) != len(s2): return 99
        return sum(a != b for a, b in zip(s1, s2))
        
    def validate_read(self, seq_or_bcs):
        # Implementation for checking extracted barcodes
        # Since splitcode extracts to FASTQ, we need to parse the FASTQ output
        # But wait, splitcode output format depends on the extraction string.
        # User config: @extract {{bc3}}<umi[10]>
        # This seems to only extract bc3 and UMI?
        # Wait, let's look at the config again.
        # It defines linker1 -> next {{bc3}}
        # linker2 -> next {{bc2}}
        # bc1_x -> next {linker2}
        # bc2_x -> next {linker1}
        # bc3_x -> next {linker1} -> wait, bc3_x next is {linker1}?
        
        # In the user config:
        # bc3_x -> next {linker1} ?? No, usually bc3 is the last one?
        # Let's re-read the structure.
        # Structure: [UMI:10][BC3:8][Linker1:30][BC2:8][Linker2:22][BC1:8][rest]
        # (Based on run_paper_benchmarks.py SplitSeqSingleEndValidityAnalyzer comments)
        
        # User config:
        # linker1 ... NEXT {{bc3}}
        # linker2 ... NEXT {{bc2}}
        # bc1_x ... NEXT {linker2}
        # bc2_x ... NEXT {linker1}
        # bc3_x ... NEXT {linker1} (Wait, this looks cyclic or weird if bc3 is at the end?)
        
        # Actually, splitcode usually parses from 5' to 3' or 3' to 5'?
        # If it's single end long read, it's usually:
        # [Adapter][BC1][Linker][BC2][Linker][BC3][UMI][PolyT]...
        # Wait, run_paper_benchmarks.py says:
        # [UMI:10][BC3:8][Linker1:30][BC2:8][Linker2:22][BC1:8]
        # So it's: 5' -> UMI -> BC3 -> L1 -> BC2 -> L2 -> BC1 -> 3'
        
        # So splitcode should find:
        # BC3 (at start?) No, UMI is first.
        # If splitcode searches:
        # The user config defines "bc3_x" entries.
        # And "@extract {{bc3}}<umi[10]>"
        
        # Let's just run it and see what it outputs.
        # We'll validate by checking if the read IDs are in the "valid set" from the raw analysis
        # OR by checking if the extracted sequence matches valid barcodes.
        pass

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

def count_fastq_reads(filepath: str) -> int:
    if not os.path.exists(filepath):
        return 0
    with open(filepath, 'r') as f:
        lines = sum(1 for _ in f)
    return lines // 4

def main():
    print(f"Testing Splitcode with NEW config on Full Dataset")
    print(f"Config: {NEW_CONFIG}")
    print("=" * 60)

    # Output file
    out_fq = "splitcode_new_out.fq"
    mapping = "splitcode_new_mapping.txt"
    
    # Clean up previous run
    if os.path.exists(out_fq): os.remove(out_fq)
    if os.path.exists(mapping): os.remove(mapping)

    # Command
    # -N 1 for single end? The config defines bc1, bc2, bc3 groups.
    # The extract string is `@extract {{bc3}}<umi[10]>`
    # This extract string seems to imply it only extracts BC3 and UMI?
    # That might be insufficient for full recovery, but let's see if it runs.
    # Also, the config provided by user has `bc1_x` ... `bc2_x` ... `bc3_x`.
    
    # We'll use -N 1 for single file input.
    cmd = f"{SPLITCODE_BIN} -c {NEW_CONFIG} --assign -N 1 -t 8 -m {mapping} -o {out_fq} {FULL_FQ}"
    
    print("Running Splitcode...")
    rt, mem, rc, stderr = run_with_memory(cmd, PROJECT_ROOT)
    
    print(f"Runtime: {rt:.2f}s")
    print(f"Memory: {mem:.2f}MB")
    print(f"Return Code: {rc}")
    
    if rc != 0:
        print("Error output:")
        print(stderr)
        
    reads = count_fastq_reads(out_fq)
    print(f"Reads extracted: {reads:,}")
    
    # Check if we got any reads
    if reads > 0:
        print("\nValidating first 5 reads...")
        with open(out_fq) as f:
            for _ in range(5):
                print(f.readline().strip()) # Header
                print(f.readline().strip()) # Seq
                f.readline()
                f.readline()
                
    # We should also check the mapping file to see what was identified
    if os.path.exists(mapping):
        print(f"\nMapping file size: {os.path.getsize(mapping)} bytes")
        with open(mapping) as f:
            valid_mapped_count = 0
            for i, line in enumerate(f):
                if i < 5: print(f"Map sample: {line.strip()}")
                valid_mapped_count += 1
            print(f"Total mapped reads in mapping file: {valid_mapped_count:,}")

if __name__ == "__main__":
    main()
