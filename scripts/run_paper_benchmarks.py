#!/usr/bin/env python3
"""
Unified Paper Benchmark Script
Runs all benchmarks for seqproc paper and outputs figures to a single folder.

Benchmarks:
1. SPLiT-seq Short-Read (Paired-End) - Raw extraction
2. SPLiT-seq Short-Read (Paired-End) - Barcode replacement (seqproc vs splitcode)
3. SPLiT-seq Long-Read (Single-End) - Raw extraction
4. 10x GridION Long-Read - Raw extraction
5. 10x PromethION Long-Read - Raw extraction

Usage:
    python scripts/run_paper_benchmarks.py --threads 4 --replicates 3
"""

import subprocess
import time
import os
import json
import tempfile
import argparse
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

# ============================================================================
# Configuration
# ============================================================================

PROJECT_ROOT = Path(os.environ.get("SEQPROC_PROJECT_ROOT", Path(__file__).parent.parent))
RESULTS_DIR = PROJECT_ROOT / "results" / "paper_figures"

# Tool binaries
SEQPROC_BIN = os.environ.get("SEQPROC_BIN", str(PROJECT_ROOT.parent / "combine-lab/seqproc/target/release/seqproc"))
MATCHBOX_BIN = os.environ.get("MATCHBOX_BIN", str(PROJECT_ROOT.parent / "matchbox/target/release/matchbox"))
SPLITCODE_BIN = os.environ.get("SPLITCODE_BIN", str(PROJECT_ROOT.parent / "splitcode/build/src/splitcode"))

# Colors
COLORS = {
    'seqproc': '#2E86AB',
    'matchbox': '#E94F37',
    'splitcode': '#7B2D8E',
}

# Dataset configurations -- file paths are overridden at startup via
# _apply_reads_level() when --reads is specified.
reads = 1000000

DATASETS = {
    # SPLiT-seq PE -- filter + edit distance (Table 2 primary config)
    'splitseq_pe_raw': {
        'name': 'SPLiT-seq PE',
        'short_name': 'SPLiT-seq PE',
        'category': 'raw',
        'r1': PROJECT_ROOT / 'data/SRR6750041_10M_R1.fastq',
        'r2': PROJECT_ROOT / 'data/SRR6750041_10M_R2.fastq',
        'mode': 'paired',
        'seqproc_geom': PROJECT_ROOT / 'configs/seqproc/splitseq_filter_edit.geom',
        'matchbox_config': PROJECT_ROOT / 'configs/matchbox/splitseq_replacement.mb',
        'splitcode_config': PROJECT_ROOT / 'configs/splitcode/splitseq_paper.config',
        'seqproc_maps': [
            PROJECT_ROOT / 'configs/seqproc/splitseq_bc3_seq2seq.tsv',
            PROJECT_ROOT / 'configs/seqproc/splitseq_bc2_seq2seq.tsv',
            PROJECT_ROOT / 'configs/seqproc/splitseq_bc1_seq2seq.tsv',
        ],
        'reads': 10_000_000,
        'tools': ['seqproc', 'matchbox', 'splitcode'],
    },

    # SPLiT-seq PE -- replacement + edit distance (supplementary config)
    'splitseq_pe_replacement': {
        'name': 'SPLiT-seq PE (Replacement)',
        'short_name': 'SPLiT-seq PE Replace',
        'category': 'replacement',
        'r1': PROJECT_ROOT / 'data/SRR6750041_10M_R1.fastq',
        'r2': PROJECT_ROOT / 'data/SRR6750041_10M_R2.fastq',
        'mode': 'paired',
        'seqproc_geom': PROJECT_ROOT / 'configs/seqproc/splitseq_replacement_edit.geom',
        'seqproc_maps': [
            PROJECT_ROOT / 'configs/seqproc/splitseq_bc3_seq2seq.tsv',
            PROJECT_ROOT / 'configs/seqproc/splitseq_bc2_seq2seq.tsv',
            PROJECT_ROOT / 'configs/seqproc/splitseq_bc1_seq2seq.tsv',
        ],
        'matchbox_config': PROJECT_ROOT / 'configs/matchbox/splitseq_replacement.mb',
        'splitcode_config': PROJECT_ROOT / 'configs/splitcode/splitseq_paper.config',
        'reads': 10_000_000,
        'tools': ['seqproc', 'matchbox', 'splitcode'],
    },
    
    # LR-SPLiT-seq (PacBio Sequel II long-read)
    'splitseq_se_raw': {
        'name': 'LR-SPLiT-seq',
        'short_name': 'LR-SPLiT-seq',
        'category': 'raw',
        'r1': PROJECT_ROOT / 'data/SRR13948564_1M.fastq',
        'r2': None,
        'mode': 'single',
        'seqproc_geom': PROJECT_ROOT / 'configs/seqproc/splitseq_singleend_primer_edit.geom',
        'matchbox_config': PROJECT_ROOT / 'configs/matchbox/splitseq_singleend.mb',
        'splitcode_config': PROJECT_ROOT / 'configs/splitcode/splitseq_singleend.config',
        'seqproc_maps': [
            PROJECT_ROOT / 'configs/seqproc/splitseq_bc3_seq2seq.tsv',
            PROJECT_ROOT / 'configs/seqproc/splitseq_bc2_seq2seq.tsv',
            PROJECT_ROOT / 'configs/seqproc/splitseq_bc1_seq2seq.tsv',
        ],
        'reads': reads,
        'tools': ['seqproc', 'matchbox', 'splitcode'],
    },
    
    '10x_short': {
        'name': '10x Chromium v2 Short Read',
        'short_name': '10x Short',
        'category': 'raw',
        'r1': PROJECT_ROOT / 'data/10x_short/SRR8315379_1M_R1.fastq',
        'r2': PROJECT_ROOT / 'data/10x_short/SRR8315379_1M_R2.fastq',
        'mode': 'paired', # While geom only uses R1, we provide both for consistency
        'seqproc_geom': PROJECT_ROOT / 'configs/seqproc/10x_v2.geom',
        'matchbox_config': PROJECT_ROOT / 'configs/matchbox/10x_v2.mb',
        'splitcode_config': PROJECT_ROOT / 'configs/splitcode/10x_v2.config',
        'reads': reads,
        'tools': ['seqproc', 'matchbox', 'splitcode'],
    },
    
    # sci-RNA-seq3
    'sciseq': {
        'name': 'sci-RNA-seq3',
        'short_name': 'sci-RNA-seq3',
        'category': 'raw',
        'r1': PROJECT_ROOT / 'data/SRR7827254_1M_1.fastq',
        'r2': PROJECT_ROOT / 'data/SRR7827254_1M_2.fastq',
        'mode': 'paired',
        'seqproc_geom': PROJECT_ROOT / 'configs/seqproc/sciseq3_edit.geom',
        'matchbox_config': PROJECT_ROOT / 'configs/matchbox/sciseq3.mb',
        'splitcode_config': PROJECT_ROOT / 'configs/splitcode/sciseq3.config',
        'reads': reads,
        'tools': ['seqproc', 'matchbox', 'splitcode'],
    },
    
    # 10x GridION Long-Read
    '10x_gridion': {
        'name': '10x GridION Long-Read',
        'short_name': '10x GridION',
        'category': 'raw',
        'r1': PROJECT_ROOT / 'data/10x/ERR9958134_1M.fastq',
        'r2': None,
        'mode': 'single',
        'seqproc_geom': PROJECT_ROOT / 'configs/seqproc/10x_longread_fwd_edit.geom',
        'seqproc_geom_rev': PROJECT_ROOT / 'configs/seqproc/10x_longread_rev_edit.geom',
        'matchbox_config': PROJECT_ROOT / 'configs/matchbox/10x_longread.mb',
        'reads': reads,
        'tools': ['seqproc', 'matchbox'],
    },
    
    # 10x PromethION Long-Read
    '10x_promethion': {
        'name': '10x PromethION Long-Read',
        'short_name': '10x PromethION',
        'category': 'raw',
        'r1': PROJECT_ROOT / 'data/10x/ERR9958135_1M.fastq',
        'r2': None,
        'mode': 'single',
        'seqproc_geom': PROJECT_ROOT / 'configs/seqproc/10x_longread_fwd_edit.geom',
        'seqproc_geom_rev': PROJECT_ROOT / 'configs/seqproc/10x_longread_rev_edit.geom',
        'matchbox_config': PROJECT_ROOT / 'configs/matchbox/10x_longread.mb',
        'reads': reads,
        'tools': ['seqproc', 'matchbox'],
    },
}


@dataclass
class BenchmarkResult:
    dataset: str
    tool: str
    runtime: float
    memory_mb: float
    reads_out: int
    reads_valid: int  # New field for valid reads
    replicate: int


# ============================================================================
# Validity Analyzer
# ============================================================================

class SplitSeqValidityAnalyzer:
    """Analyzes SPLiT-seq reads for validity (d<=1)."""
    
    LINKER1 = "GTGGCCGCTGTTTCGCATCGGCGTACGACT"  # 30bp
    LINKER2 = "ATCCACGTGCTTGAGAGGCCAGAGCATTCG"  # 30bp
    
    def __init__(self, bc1_map, bc2_map, bc3_map):
        self.bc1_wl = self._load_whitelist(bc1_map, truncate=6)
        self.bc2_wl = self._load_whitelist(bc2_map)
        self.bc3_wl = self._load_whitelist(bc3_map)
        self.valid_ids = set()
        
    def _load_whitelist(self, path, truncate=None):
        wl = set()
        with open(path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    seq = parts[1]
                    if truncate and len(seq) > truncate:
                        seq = seq[:truncate]
                    wl.add(seq)
        return wl
        
    def _hamming(self, s1, s2):
        if len(s1) != len(s2): return 99
        return sum(a != b for a, b in zip(s1, s2))
    
    def _check_wl(self, bc, wl):
        """Check if barcode matches whitelist within Hamming distance 1."""
        if bc in wl: return True
        for cand in wl:
            if self._hamming(bc, cand) <= 1: return True
        return False
        
    def _find_linker(self, read, linker, start=0):
        best_pos, best_dist = -1, 100
        # Optimization: Only search in plausible window
        # L1 usually around 18, L2 around 56
        search_end = min(len(read) - len(linker) + 1, start + 40)
        
        for i in range(start, search_end):
            dist = self._hamming(read[i:i+len(linker)], linker)
            if dist < best_dist:
                best_dist = dist
                best_pos = i
                if dist <= 1: break # optimization
        return best_pos, best_dist
        
    def analyze_fastqs(self, r1_path, r2_path):
        """Analyze raw FASTQs and return set of valid read IDs (d<=1)."""
        print(f"    Analyzing raw input for validity (this may take 15-20s)...")
        valid_ids = set()
        
        # We only need R2 for validity
        with open(r2_path, 'r') as f:
            while True:
                header = f.readline()
                if not header: break
                seq = f.readline().strip()
                f.readline()
                f.readline()
                
                read_id = header.strip().split()[0].replace('@', '')
                
                # Check structure
                l1_pos, l1_dist = self._find_linker(seq, self.LINKER1, 0)
                if l1_dist > 3: continue
                
                l2_pos, l2_dist = self._find_linker(seq, self.LINKER2, l1_pos + 30)
                if l2_dist > 3: continue
                
                # Extract BCs
                # NN(2) + UMI(10) + BC3(8) + L1(30) + BC2(8) + L2(30) + BC1(6)
                bc3 = seq[l1_pos-8:l1_pos]
                bc2 = seq[l1_pos+30:l1_pos+38]
                bc1 = seq[l2_pos+30:l2_pos+36]
                
                if len(bc3) != 8 or len(bc2) != 8 or len(bc1) != 6:
                    continue
                    
                # Check validity (d<=1)
                # Optimization: check exact match first (fastest)
                if bc3 in self.bc3_wl and bc2 in self.bc2_wl and bc1 in self.bc1_wl:
                    valid_ids.add(read_id)
                    continue
                    
                # Check d<=1
                # This is slow in Python, but we only do it for non-exact matches
                if self._check_wl(bc3, self.bc3_wl) and self._check_wl(bc2, self.bc2_wl) and self._check_wl(bc1, self.bc1_wl):
                    valid_ids.add(read_id)
                    
        print(f"    Found {len(valid_ids):,} valid reads (d<=1) in input.")
        self.valid_ids = valid_ids
        return valid_ids


class SplitSeqSingleEndValidityAnalyzer:
    """Analyzes SPLiT-seq Single-End reads for validity (d<=1)."""
    
    # SE Linkers (different from PE)
    LINKER1 = "GTGGCCGATGTTTCGCATCGGCGTACGACT"  # 30bp
    # Linker 2 is skipped in geometry but present in sequence: ATCCACGTGCTTGAGACTGTGG (22bp)
    LINKER2 = "ATCCACGTGCTTGAGACTGTGG" 
    
    def __init__(self, bc1_map, bc2_map, bc3_map):
        self.bc1_wl = self._load_whitelist(bc1_map)
        self.bc2_wl = self._load_whitelist(bc2_map)
        self.bc3_wl = self._load_whitelist(bc3_map)
        self.valid_ids = set()
        
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
    
    def _check_wl(self, bc, wl):
        """Check if barcode matches whitelist within Hamming distance 1."""
        if bc in wl: return True
        for cand in wl:
            if self._hamming(bc, cand) <= 1: return True
        return False
        
    def _find_linker(self, read, linker, start=0):
        best_pos, best_dist = -1, 100
        search_end = min(len(read) - len(linker) + 1, start + 40)
        
        for i in range(start, search_end):
            dist = self._hamming(read[i:i+len(linker)], linker)
            if dist < best_dist:
                best_dist = dist
                best_pos = i
                if dist <= 1: break
        return best_pos, best_dist
        
    def analyze_fastqs(self, r1_path):
        """Analyze raw R1 FASTQ and return set of valid read IDs (d<=1)."""
        print(f"    Analyzing raw input for validity (Single-End)...")
        valid_ids = set()
        
        with open(r1_path, 'r') as f:
            while True:
                header = f.readline()
                if not header: break
                seq = f.readline().strip()
                f.readline()
                f.readline()
                
                read_id = header.strip().split()[0].replace('@', '')
                
                # Structure: [UMI:10][BC3:8][Linker1:30][BC2:8][Linker2:22][BC1:8]
                # Find Linker 1
                l1_pos, l1_dist = self._find_linker(seq, self.LINKER1, 10) # Start search after UMI+BC3 approx
                if l1_dist > 3: continue
                
                # Check bounds
                if l1_pos < 18: continue # Needs 10+8 before it
                
                bc3 = seq[l1_pos-8:l1_pos]
                bc2 = seq[l1_pos+30:l1_pos+38]
                bc1 = seq[l1_pos+30+8+22:l1_pos+30+8+22+8] # L1(30) + BC2(8) + L2(22)
                
                if len(bc3) != 8 or len(bc2) != 8 or len(bc1) != 8:
                    continue

                if self._check_wl(bc3, self.bc3_wl) and self._check_wl(bc2, self.bc2_wl) and self._check_wl(bc1, self.bc1_wl):
                    valid_ids.add(read_id)
                    
        print(f"    Found {len(valid_ids):,} valid reads (d<=1) in input.")
        self.valid_ids = valid_ids
        return valid_ids


class SciSeqValidityAnalyzer:
    """Analyzes Sci-Seq reads for validity (Anchor check)."""
    
    ANCHOR = "CAGAGC"
    
    def __init__(self):
        self.valid_ids = set()
        
    def _hamming(self, s1, s2):
        if len(s1) != len(s2): return 99
        return sum(a != b for a, b in zip(s1, s2))
        
    def analyze_fastqs(self, r1_path, r2_path=None):
        """Analyze R1 FASTQ for validity."""
        print(f"    Analyzing Sci-Seq input for validity (Anchor Check)...")
        valid_ids = set()
        
        with open(r1_path, 'r') as f:
            while True:
                header = f.readline()
                if not header: break
                seq = f.readline().strip()
                f.readline()
                f.readline()
                
                # Check for anchor CAGAGC
                # BC1 is 9-10bp, so Anchor should be around index 9 or 10
                # Search window: 9 to 11 (0-indexed)
                found = False
                for i in range(8, 12): 
                    if i + len(self.ANCHOR) <= len(seq):
                        dist = self._hamming(seq[i:i+len(self.ANCHOR)], self.ANCHOR)
                        if dist <= 1:
                            found = True
                            break
                
                if found:
                    valid_ids.add(header.strip().split()[0].replace('@', ''))
                    
        print(f"    Found {len(valid_ids):,} valid reads in input.")
        self.valid_ids = valid_ids
        return valid_ids


class TenXValidityAnalyzer:
    """Analyzes 10x reads for validity (Structure check)."""
    
    # 10x v2/v3 Primer (for Long Reads)
    PRIMER = "CTACACGACGCTCTTCCGATCT"
    
    def __init__(self, is_short_read=False):
        self.valid_ids = set()
        self.is_short_read = is_short_read
        self.whitelist = self._load_whitelist()
        
    def _load_whitelist(self):
        # Prefer v3 whitelist if available
        wl_paths = [
            "/home/ubuntu/3M-february-2018.txt.gz",
            "/home/ubuntu/737K-august-2016.txt"
        ]
        wl = set()
        for path in wl_paths:
            if os.path.exists(path):
                print(f"    Loading 10x whitelist from {path}...")
                import gzip
                opener = gzip.open if path.endswith('.gz') else open
                mode = 'rt' if path.endswith('.gz') else 'r'
                with opener(path, mode) as f:
                    for line in f:
                        bc = line.strip().split(',')[0] # handle CSV if needed
                        if len(bc) >= 16:
                            wl.add(bc)
                print(f"    Loaded {len(wl):,} barcodes.")
                break
        return wl
        
    def _hamming(self, s1, s2):
        if len(s1) != len(s2): return 99
        return sum(a != b for a, b in zip(s1, s2))
    
    def _rc(self, seq):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        return "".join(complement.get(base, base) for base in reversed(seq))

    def _find_primer_and_barcode(self, seq, primer):
        # Search for primer in forward
        primer_len = len(primer)
        best_dist = 99
        
        # Heuristic search window: usually at start or end depending on strand
        search_regions = [(0, 100), (max(0, len(seq)-100), len(seq))]
        
        # Forward Search: Primer -> Barcode
        # Structure: ...[Primer][Barcode]...
        for start, end in search_regions:
            for i in range(start, min(end, len(seq) - primer_len - 16 + 1)):
                dist = self._hamming(seq[i:i+primer_len], primer)
                if dist <= 3:
                    # Found primer, extract barcode after
                    bc = seq[i+primer_len : i+primer_len+16]
                    if len(bc) == 16:
                        return bc
        
        # Reverse Complement Search
        # Structure: ...[Barcode_RC][Primer_RC]...
        # Real barcode is RC(Barcode_RC)
        primer_rc = self._rc(primer)
        primer_rc_len = len(primer_rc)
        
        for start, end in search_regions:
            # Need 16bp before primer
            search_start = max(start, 16)
            for i in range(search_start, min(end, len(seq) - primer_rc_len + 1)):
                dist = self._hamming(seq[i:i+primer_rc_len], primer_rc)
                if dist <= 3:
                    # Found primer_rc, extract barcode before
                    bc_rc = seq[i-16 : i]
                    if len(bc_rc) == 16:
                        return self._rc(bc_rc)
                
        return None
        
    def analyze_fastqs(self, r1_path, r2_path=None):
        """Analyze FASTQ for validity."""
        mode = "Short Read (Length=26)" if self.is_short_read else "Long Read (Primer Check + Whitelist)"
        print(f"    Analyzing 10x input for validity ({mode})...")
        valid_ids = set()
        
        # Only analyze if we are in short read mode (structure) or have whitelist
        if not self.is_short_read and not self.whitelist:
             print("    [WARNING] No whitelist found for 10x validation. Skipping validity check.")
             return valid_ids

        with open(r1_path, 'r') as f:
            while True:
                header = f.readline()
                if not header: break
                seq = f.readline().strip()
                f.readline()
                f.readline()
                
                read_id = header.strip().split()[0].replace('@', '')

                if self.is_short_read:
                    # 10x v2 R1 must be exactly 26bp (16 BC + 10 UMI)
                    if len(seq) == 26:
                        valid_ids.add(read_id)
                else:
                    # Long Read: Find primer, extract BC, check whitelist
                    bc = self._find_primer_and_barcode(seq, self.PRIMER)
                    if bc and bc in self.whitelist:
                        valid_ids.add(read_id)
                    
        print(f"    Found {len(valid_ids):,} valid reads in input (Gold Standard).")
        self.valid_ids = valid_ids
        return valid_ids

# ============================================================================
# Helper Functions
# ============================================================================

def run_with_memory(cmd: str, cwd=None) -> Tuple[float, float, int, str]:
    """Run command and return (runtime, peak_memory_mb, returncode, stderr)."""
    time_cmd = f"/usr/bin/time -v {cmd}"
    start = time.time()
    result = subprocess.run(time_cmd, shell=True, capture_output=True, text=True, cwd=cwd)
    runtime = time.time() - start
    
    peak_mem_kb = 0
    for line in result.stderr.split('\n'):
        if 'Maximum resident set size' in line:
            peak_mem_kb = int(line.split(':')[1].strip())
            break
    
    return runtime, peak_mem_kb / 1024, result.returncode, result.stderr


def count_fastq_reads(filepath: str) -> int:
    """Count reads in FASTQ file."""
    if not os.path.exists(filepath):
        return 0
    with open(filepath, 'r') as f:
        lines = sum(1 for _ in f)
    return lines // 4


def count_tsv_lines(filepath: str) -> int:
    """Count lines in TSV file."""
    if not os.path.exists(filepath):
        return 0
    with open(filepath, 'r') as f:
        return sum(1 for _ in f)


# ============================================================================
# Tool Runners
# ============================================================================

def run_seqproc(dataset: dict, tmpdir: str, threads: int) -> Tuple[float, float, int]:
    """Run seqproc on dataset."""
    out1 = f"{tmpdir}/seqproc_R1.fq"
    out2 = f"{tmpdir}/seqproc_R2.fq"
    
    # Check if we need dual geometry (for 10x long read mixed orientation)
    if 'seqproc_geom_rev' in dataset:
        out1_fwd = f"{tmpdir}/seqproc_R1_fwd.fq"
        out1_rev = f"{tmpdir}/seqproc_R1_rev.fq"
        
        # Run forward pass (single output)
        cmd_fwd = f"{SEQPROC_BIN} --geom {dataset['seqproc_geom']} --file1 {dataset['r1']} --out1 {out1_fwd} --threads {threads}"
        rt1, mem1, rc1, stderr1 = run_with_memory(cmd_fwd, PROJECT_ROOT)
        if rc1 != 0:
            print(f"\n[ERROR] Seqproc fwd pass failed with rc={rc1}:\n{stderr1}")
        
        # Run reverse pass (single output)
        cmd_rev = f"{SEQPROC_BIN} --geom {dataset['seqproc_geom_rev']} --file1 {dataset['r1']} --out1 {out1_rev} --threads {threads}"
        rt2, mem2, rc2, stderr2 = run_with_memory(cmd_rev, PROJECT_ROOT)
        if rc2 != 0:
            print(f"\n[ERROR] Seqproc rev pass failed with rc={rc2}:\n{stderr2}")
        
        # Merge outputs
        if os.path.exists(out1_fwd) or os.path.exists(out1_rev):
            with open(out1, 'wb') as outfile:
                if os.path.exists(out1_fwd):
                    with open(out1_fwd, 'rb') as f: outfile.write(f.read())
                if os.path.exists(out1_rev):
                    with open(out1_rev, 'rb') as f: outfile.write(f.read())
            
        runtime = rt1 + rt2
        memory_mb = max(mem1, mem2)
        
        reads = count_fastq_reads(out1)
        return runtime, memory_mb, reads

    if dataset['mode'] == 'single':
        cmd = f"{SEQPROC_BIN} --geom {dataset['seqproc_geom']} --file1 {dataset['r1']} --out1 {out1} --threads {threads}"
    else:
        cmd = f"{SEQPROC_BIN} --geom {dataset['seqproc_geom']} --file1 {dataset['r1']} --file2 {dataset['r2']} --out1 {out1} --out2 {out2} --threads {threads}"
    
    # Add map files for replacement benchmark
    if 'seqproc_maps' in dataset:
        for map_file in dataset['seqproc_maps']:
            cmd += f" -a {map_file}"
    
    runtime, memory_mb, rc, stderr = run_with_memory(cmd, PROJECT_ROOT)
    
    if rc != 0:
        print(f"\n[ERROR] Seqproc failed with rc={rc}:\n{stderr}")
    
    output_file = out2 if dataset['mode'] == 'paired' else out1
    reads = count_fastq_reads(output_file)
    
    return runtime, memory_mb, reads


def run_matchbox(dataset: dict, tmpdir: str, threads: int) -> Tuple[float, float, int]:
    """Run matchbox on dataset."""
    out_tsv = f"{tmpdir}/matchbox_out.tsv"
    
    # Clean stale matchbox output files from previous runs
    for fq in ['mb_r1.fq', 'mb_r2.fq']:
        stale = PROJECT_ROOT / fq
        if stale.exists():
            stale.unlink()
    
    # Determine input arguments and output expectations
    args = ""
    if dataset['name'] == '10x Chromium v2 Short Read':
         # 10x Short is Paired but we only scan R1 in the raw benchmark
         input_file = dataset['r1']
         args = str(input_file)
    elif dataset['mode'] == 'paired':
         # Pass both files for paired datasets (SciSeq, SPLiT-seq) using -p flag
         args = f"{dataset['r1']} -p {dataset['r2']}"
    else:
         # Single end
         args = str(dataset['r1'])
        
    cmd = f'{MATCHBOX_BIN} -e 0.2 -t {threads} -s {dataset["matchbox_config"]} {args} > {out_tsv}'
    
    # Run matchbox from PROJECT_ROOT so relative config paths work
    runtime, memory_mb, rc, stderr = run_with_memory(cmd, PROJECT_ROOT)
    
    if rc != 0:
        print(f"\n[ERROR] Matchbox failed with rc={rc}:\n{stderr}")
    
    # Check for generated FASTQ files (mb_r1.fq, mb_r2.fq) in PROJECT_ROOT
    # and move them to tmpdir
    reads = 0
    generated_fastqs = []
    
    for fq in ['mb_r1.fq', 'mb_r2.fq']:
        src = PROJECT_ROOT / fq
        # Remove stale files from previous runs before checking
        dst = Path(tmpdir) / fq
        if src.exists():
            os.rename(src, dst)
            generated_fastqs.append(dst)
            
    if generated_fastqs:
        # Count reads from the first generated FASTQ
        reads = count_fastq_reads(str(generated_fastqs[0]))
    else:
        # Fallback to TSV line count
        reads = count_tsv_lines(out_tsv)
        
    if reads == 0 and rc == 0:
         # Only warn if it's NOT the raw extraction benchmark which outputs TSV to stdout
         # The raw benchmarks redirect stdout to out_tsv, so reads count should come from out_tsv
         if not generated_fastqs and count_tsv_lines(out_tsv) == 0:
            print(f"\n[WARNING] Matchbox produced 0 reads. Stderr:\n{stderr}")
    
    return runtime, memory_mb, reads


def run_splitcode(dataset: dict, tmpdir: str, threads: int) -> Tuple[float, float, int]:
    """Run splitcode on dataset."""
    mapping = f"{tmpdir}/splitcode_mapping.txt"
    
    # Detect whether config has tag lines (needs --assign) or only @extract directives
    with open(dataset['splitcode_config']) as f:
        config_text = f.read()
    has_tags = False
    for line in config_text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith('#') or stripped.startswith('@'):
            continue
        # Header row (ID/group etc.) or actual tag data rows are tab-separated
        fields = stripped.split('\t')
        if len(fields) >= 3 and fields[0] not in ('ID', 'group'):
            has_tags = True
            break
    assign_flag = "--assign" if has_tags else ""
    mapping_flag = f"-m {mapping}" if has_tags else ""
    
    if dataset['mode'] == 'single':
        out_fq = f"{tmpdir}/splitcode_out.fq"
        cmd = f"{SPLITCODE_BIN} -c {dataset['splitcode_config']} {assign_flag} -N 1 -t {threads} {mapping_flag} -o {out_fq} {dataset['r1']}"
    else:
        out1 = f"{tmpdir}/splitcode_R1.fq"
        out2 = f"{tmpdir}/splitcode_R2.fq"
        cmd = f"{SPLITCODE_BIN} -c {dataset['splitcode_config']} {assign_flag} -N 2 -t {threads} {mapping_flag} -o {out1},{out2} {dataset['r1']} {dataset['r2']}"
    
    runtime, memory_mb, rc, stderr = run_with_memory(cmd, PROJECT_ROOT)
    
    if rc != 0:
        print(f"\n[ERROR] Splitcode failed with rc={rc}:\n{stderr}")
    
    if dataset['mode'] == 'single':
        reads = count_fastq_reads(out_fq)
    else:
        reads = count_fastq_reads(out2)
        
    if reads == 0 and rc == 0:
        print(f"\n[WARNING] Splitcode produced 0 reads. Stderr:\n{stderr}")
    
    return runtime, memory_mb, reads


# ============================================================================
# Benchmark Runner
# ============================================================================

def run_benchmarks(threads: int, replicates: int, dataset_filter=None) -> List[BenchmarkResult]:
    """Run all benchmarks."""
    results = []
    
    for dataset_key, dataset in DATASETS.items():
        if dataset_filter and dataset_key not in dataset_filter:
            continue
        # Check if data exists
        if not os.path.exists(dataset['r1']):
            print(f"  [SKIP] Data not found: {dataset['r1']}")
            continue
        
        print(f"\n{'='*60}")
        print(f"Dataset: {dataset['name']}")
        print(f"{'='*60}")
        
        # Instantiate validity analyzer ONCE per dataset (outside replicate loop)
        validity_analyzer = None
        
        if 'splitseq' in dataset_key:
            # SPLiT-seq datasets
            if 'seqproc_maps' in dataset:
                if dataset['mode'] == 'paired':
                    # PE Raw and PE Replacement
                    validity_analyzer = SplitSeqValidityAnalyzer(
                        dataset['seqproc_maps'][2],  # bc1
                        dataset['seqproc_maps'][1],  # bc2
                        dataset['seqproc_maps'][0]   # bc3
                    )
                    validity_analyzer.analyze_fastqs(str(dataset['r1']), str(dataset['r2']))
                else:
                    # SE Long Read
                    validity_analyzer = SplitSeqSingleEndValidityAnalyzer(
                        dataset['seqproc_maps'][2],  # bc1
                        dataset['seqproc_maps'][1],  # bc2
                        dataset['seqproc_maps'][0]   # bc3
                    )
                    validity_analyzer.analyze_fastqs(str(dataset['r1']))
        
        elif dataset_key.startswith('10x_'):
            # 10x datasets
            is_short = (dataset_key == '10x_short')
            validity_analyzer = TenXValidityAnalyzer(is_short_read=is_short)
            validity_analyzer.analyze_fastqs(str(dataset['r1']))
        
        elif dataset_key == 'sciseq':
            # Sci-Seq dataset
            validity_analyzer = SciSeqValidityAnalyzer()
            validity_analyzer.analyze_fastqs(str(dataset['r1']))
        
        for rep in range(1, replicates + 1):
            print(f"\n  Replicate {rep}/{replicates}:")
            
            with tempfile.TemporaryDirectory() as tmpdir:
                for tool in dataset['tools']:
                    print(f"    {tool}...", end=" ", flush=True)
                    
                    if tool == 'seqproc':
                        runtime, memory, reads = run_seqproc(dataset, tmpdir, threads)
                    elif tool == 'matchbox':
                        runtime, memory, reads = run_matchbox(dataset, tmpdir, threads)
                    elif tool == 'splitcode':
                        runtime, memory, reads = run_splitcode(dataset, tmpdir, threads)
                    else:
                        continue
                    
                    # Calculate valid reads if applicable
                    reads_valid = 0
                    if validity_analyzer:
                        output_ids = set()
                        # Get output read IDs
                        if tool == 'seqproc':
                            # Handle both paired (R2) and single (R1) output
                            if dataset['mode'] == 'paired':
                                out_file = f"{tmpdir}/seqproc_R2.fq"
                            else:
                                out_file = f"{tmpdir}/seqproc_R1.fq"
                                
                            if os.path.exists(out_file):
                                with open(out_file) as f:
                                    while True:
                                        h = f.readline()
                                        if not h: break
                                        f.readline(); f.readline(); f.readline()
                                        output_ids.add(h.strip().split()[0].replace('@', ''))

                        elif tool == 'splitcode':
                            if dataset['mode'] == 'paired':
                                out_file = f"{tmpdir}/splitcode_R2.fq"
                            else:
                                out_file = f"{tmpdir}/splitcode_out.fq"
                                
                            if os.path.exists(out_file):
                                with open(out_file) as f:
                                    while True:
                                        h = f.readline()
                                        if not h: break
                                        f.readline(); f.readline(); f.readline()
                                        output_ids.add(h.strip().split()[0].replace('@', ''))

                        elif tool == 'matchbox':
                            mb_fq_files = [f for f in [f"{tmpdir}/mb_r1.fq", f"{tmpdir}/mb_r2.fq"] if os.path.exists(f)]
                            
                            if mb_fq_files:
                                with open(mb_fq_files[0]) as f:
                                    while True:
                                        h = f.readline()
                                        if not h: break
                                        f.readline(); f.readline(); f.readline()
                                        output_ids.add(h.strip().split()[0].replace('@', ''))
                            else:
                                out_file = f"{tmpdir}/matchbox_out.tsv"
                                if os.path.exists(out_file):
                                    with open(out_file) as f:
                                        for line in f:
                                            parts = line.strip().split('\t')
                                            if parts:
                                                output_ids.add(parts[0])

                        reads_valid = len(output_ids.intersection(validity_analyzer.valid_ids))
                        print(f"{runtime:.2f}s, {memory:.1f}MB, {reads:,} reads ({reads_valid:,} valid)")
                    else:
                        print(f"{runtime:.2f}s, {memory:.1f}MB, {reads:,} reads")
                        
                    results.append(BenchmarkResult(
                        dataset=dataset_key,
                        tool=tool,
                        runtime=runtime,
                        memory_mb=memory,
                        reads_out=reads,
                        reads_valid=reads_valid,
                        replicate=rep
                    ))
    
    return results


# ============================================================================
# Figure Generation
# ============================================================================

def generate_all_figures(results: List[BenchmarkResult], output_dir: Path):
    """Generate all paper figures to a single folder."""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Aggregate results
    data = defaultdict(lambda: defaultdict(list))
    for r in results:
        data[r.dataset][r.tool].append({
            'runtime': r.runtime,
            'memory': r.memory_mb,
            'reads': r.reads_out
        })
    
    # ========================================================================
    # Figure 1: Performance Distribution (Violin Plots)
    # ========================================================================
    # Combine all datasets into one figure with subplots
    all_datasets = [k for k in DATASETS.keys() if k in data]
    
    if all_datasets:
        fig, axes = plt.subplots(len(all_datasets), 2, figsize=(12, 4 * len(all_datasets)),
                                  squeeze=False)
        
        for idx, ds_key in enumerate(all_datasets):
            ds_info = DATASETS[ds_key]
            ax_runtime = axes[idx][0]
            ax_memory = axes[idx][1]
            
            tools = [t for t in ['seqproc', 'matchbox', 'splitcode'] if t in data[ds_key]]
            runtime_data = [[r['runtime'] for r in data[ds_key][t]] for t in tools]
            memory_data = [[r['memory'] for r in data[ds_key][t]] for t in tools]
            
            # Runtime Violin
            if runtime_data:
                parts = ax_runtime.violinplot(runtime_data, showmeans=True, showextrema=True)
                for i, pc in enumerate(parts['bodies']):
                    pc.set_facecolor(COLORS[tools[i]])
                    pc.set_alpha(0.7)
                
                ax_runtime.set_xticks(np.arange(1, len(tools) + 1))
                ax_runtime.set_xticklabels([t.capitalize() for t in tools])
            ax_runtime.set_ylabel('Runtime (s)')
            ax_runtime.set_title(f"{ds_info['short_name']} - Runtime")
            ax_runtime.grid(axis='y', alpha=0.3)
            
            # Memory Violin
            if memory_data:
                parts = ax_memory.violinplot(memory_data, showmeans=True, showextrema=True)
                for i, pc in enumerate(parts['bodies']):
                    pc.set_facecolor(COLORS[tools[i]])
                    pc.set_alpha(0.7)
                
                ax_memory.set_xticks(np.arange(1, len(tools) + 1))
                ax_memory.set_xticklabels([t.capitalize() for t in tools])
            ax_memory.set_ylabel('Peak Memory (MB)')
            ax_memory.set_title(f"{ds_info['short_name']} - Memory")
            ax_memory.grid(axis='y', alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_dir / 'fig1_performance_distribution.pdf', bbox_inches='tight')
        plt.savefig(output_dir / 'fig1_performance_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: fig1_performance_distribution.png")

    # ========================================================================
    # Figure 2: Read Recovery Table
    # ========================================================================
    if all_datasets:
        # Pre-build valid reads lookup to avoid O(n^2) scan
        valid_lookup = defaultdict(list)
        for r in results:
            valid_lookup[(r.dataset, r.tool)].append(r.reads_valid)
        
        # Prepare table data
        columns = ['Dataset', 'Tool', 'Total Input', 'Output Reads', 'Recovery %', 'Valid Reads (d<=1)', 'Valid %', 'Filtered Out']
        cell_text = []
        
        for ds_key in all_datasets:
            ds_info = DATASETS[ds_key]
            total_input = ds_info['reads']
            
            for tool in ['seqproc', 'matchbox', 'splitcode']:
                if tool in data[ds_key]:
                    runs = data[ds_key][tool]
                    mean_out = np.mean([r['reads'] for r in runs])
                    
                    # Get valid counts if available
                    valid_counts = valid_lookup.get((ds_key, tool), [])
                    if any(v > 0 for v in valid_counts):
                        mean_valid = np.mean(valid_counts)
                        valid_pct = (mean_valid / total_input) * 100
                        valid_str = f"{mean_valid:,.0f}"
                        valid_pct_str = f"{valid_pct:.1f}%"
                    else:
                        valid_str = "N/A"
                        valid_pct_str = "-"

                    filtered = total_input - mean_out
                    
                    cell_text.append([
                        ds_info['short_name'],
                        tool.capitalize(),
                        f"{total_input:,}",
                        f"{mean_out:,.0f}",
                        f"{(mean_out/total_input)*100:.1f}%",
                        valid_str,
                        valid_pct_str,
                        f"{filtered:,.0f}"
                    ])

        # Create table plot
        # Approx height: header + rows
        fig, ax = plt.subplots(figsize=(14, len(cell_text) * 0.4 + 1.5))
        ax.axis('off')
        
        table = ax.table(cellText=cell_text, colLabels=columns, loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.5)
        
        # Style header
        for (row, col), cell in table.get_celld().items():
            if row == 0:
                cell.set_text_props(weight='bold')
                cell.set_facecolor('#f0f0f0')
            cell.set_edgecolor('black')
        
        plt.title("Read Recovery & Validity Statistics", fontweight='bold', pad=10)
        plt.tight_layout()
        plt.savefig(output_dir / 'fig2_recovery_table.pdf', bbox_inches='tight')
        plt.savefig(output_dir / 'fig2_recovery_table.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: fig2_recovery_table.png")

    # ========================================================================
    # Figure 3: Summary Table
    # ========================================================================
    if all_datasets:
        columns = ['Dataset', 'Tool', 'Mean Runtime (s)', 'Mean Memory (MB)']
        cell_text = []
        
        for ds_key in all_datasets:
            ds_info = DATASETS[ds_key]
            for tool in ['seqproc', 'matchbox', 'splitcode']:
                if tool in data[ds_key]:
                    runs = data[ds_key][tool]
                    mean_rt = np.mean([r['runtime'] for r in runs])
                    std_rt = np.std([r['runtime'] for r in runs])
                    mean_mem = np.mean([r['memory'] for r in runs])
                    std_mem = np.std([r['memory'] for r in runs])
                    
                    cell_text.append([
                        ds_info['short_name'],
                        tool.capitalize(),
                        f"{mean_rt:.2f} ± {std_rt:.2f}",
                        f"{mean_mem:.1f} ± {std_mem:.1f}"
                    ])
        
        fig, ax = plt.subplots(figsize=(10, len(cell_text) * 0.4 + 1.5))
        ax.axis('off')
        
        table = ax.table(cellText=cell_text, colLabels=columns, loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.5)
        
        # Style header
        for (row, col), cell in table.get_celld().items():
            if row == 0:
                cell.set_text_props(weight='bold')
                cell.set_facecolor('#e6e6fa')
            cell.set_edgecolor('black')

        plt.title("Performance Summary", fontweight='bold', pad=10)
        plt.tight_layout()
        plt.savefig(output_dir / 'fig3_summary_table.pdf', bbox_inches='tight')
        plt.savefig(output_dir / 'fig3_summary_table.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: fig3_summary_table.png")


def save_results_json(results: List[BenchmarkResult], output_dir: Path):
    """Save results as JSON."""
    data = defaultdict(lambda: defaultdict(list))
    for r in results:
        data[r.dataset][r.tool].append({
            'runtime': r.runtime,
            'memory_mb': r.memory_mb,
            'reads_out': r.reads_out,
            'reads_valid': r.reads_valid,
            'replicate': r.replicate
        })
    
    summary = {}
    for ds_key, tool_data in data.items():
        summary[ds_key] = {
            'name': DATASETS[ds_key]['name'],
            'total_reads': DATASETS[ds_key]['reads'],
            'tools': {}
        }
        for tool, runs in tool_data.items():
            mean_valid = float(np.mean([r['reads_valid'] for r in runs]))
            summary[ds_key]['tools'][tool] = {
                'mean_runtime': float(np.mean([r['runtime'] for r in runs])),
                'std_runtime': float(np.std([r['runtime'] for r in runs])),
                'mean_memory_mb': float(np.mean([r['memory_mb'] for r in runs])),
                'mean_reads_out': int(np.mean([r['reads_out'] for r in runs])),
                'mean_reads_valid': int(mean_valid),
                'recovery_rate': float(np.mean([r['reads_out'] for r in runs]) / DATASETS[ds_key]['reads'] * 100),
                'valid_rate': float(mean_valid / DATASETS[ds_key]['reads'] * 100) if mean_valid > 0 else 0.0
            }
    
    with open(output_dir / 'benchmark_results.json', 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"  Saved: benchmark_results.json")


# ============================================================================
# Main
# ============================================================================

def _apply_reads_level(reads_level: str):
    """Override DATASETS file paths and read counts based on --reads level.

    Maps the benchmark DATASETS keys to data_config's canonical dataset keys
    and replaces r1, r2, reads for each matching entry.
    """
    from data_config import resolve_datasets as _resolve

    resolved = _resolve(reads_level)

    # Mapping from benchmark DATASETS keys to data_config canonical keys.
    # Multiple benchmark entries can map to the same canonical dataset
    # (e.g. splitseq_pe_raw and splitseq_pe_replacement both use SRR6750041).
    _KEY_MAP = {
        "splitseq_pe_raw": "splitseq_pe",
        "splitseq_pe_replacement": "splitseq_pe",
        "splitseq_se_raw": "lr_splitseq",
        "10x_short": "10x_short",
        "sciseq": "sciseq",
        # gridion/promethion are not in data_config (supplementary); skip them.
    }

    for bench_key, canon_key in _KEY_MAP.items():
        if bench_key in DATASETS and canon_key in resolved:
            src = resolved[canon_key]
            DATASETS[bench_key]["r1"] = src["r1"]
            DATASETS[bench_key]["r2"] = src["r2"]
            DATASETS[bench_key]["reads"] = src["reads"]


def main():
    parser = argparse.ArgumentParser(description='Run all paper benchmarks')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads')
    parser.add_argument('--replicates', type=int, default=3, help='Number of replicates')
    parser.add_argument('--output', type=str, default=None, help='Output directory')
    parser.add_argument('--datasets', type=str, nargs='+', default=None,
                        help='Specific dataset keys to run (default: all)')
    parser.add_argument('--reads', type=str, choices=['1m', 'full'], default='1m',
                        help="Dataset size: '1m' (default) or 'full'")
    args = parser.parse_args()

    # Override file paths for the requested reads level
    _apply_reads_level(args.reads)
    
    output_dir = Path(args.output) if args.output else RESULTS_DIR
    
    print("=" * 70)
    print("SEQPROC PAPER BENCHMARKS")
    print("=" * 70)
    print(f"Threads: {args.threads}")
    print(f"Replicates: {args.replicates}")
    print(f"Output: {output_dir}")
    
    if args.datasets:
        print(f"Datasets: {', '.join(args.datasets)}")
    
    # Run benchmarks
    results = run_benchmarks(args.threads, args.replicates, dataset_filter=args.datasets)
    
    if not results:
        print("\nNo benchmark results collected. Check that data files exist.")
        return
    
    # Generate figures
    print("\n" + "=" * 70)
    print("GENERATING FIGURES")
    print("=" * 70)
    generate_all_figures(results, output_dir)
    save_results_json(results, output_dir)
    
    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)
    print(f"All figures saved to: {output_dir}")


if __name__ == "__main__":
    main()
