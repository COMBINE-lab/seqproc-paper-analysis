#!/usr/bin/env python3
"""
SPLiT-seq Read Discrepancy Investigation

This script performs a deep analysis of why seqproc and splitcode produce
different read counts. It examines:
1. Raw read recovery from both tools
2. Barcode validity (within distance 1 of whitelist)
3. Per-barcode breakdown of matches/mismatches
4. Detailed comparison of accepted vs rejected reads

Goal: 100% insight into the read count discrepancy.
"""

import subprocess
import os
import sys
import tempfile
from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Set
import time

# ============================================================================
# Configuration
# ============================================================================

PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
CONFIG_DIR = PROJECT_ROOT / "configs"

# Binaries
SEQPROC_BIN = str(PROJECT_ROOT.parent / "seqproc/target/release/seqproc")
SPLITCODE_BIN = str(PROJECT_ROOT.parent / "splitcode/build/src/splitcode")

# Data files (1M subset for quick analysis)
R1_FILE = DATA_DIR / "SRR6750041_1M_R1.fastq"
R2_FILE = DATA_DIR / "SRR6750041_1M_R2.fastq"

# Whitelist files
BC1_WHITELIST = CONFIG_DIR / "seqproc/splitseq_bc1_map.tsv"
BC2_WHITELIST = CONFIG_DIR / "seqproc/splitseq_bc2_map.tsv"
BC3_WHITELIST = CONFIG_DIR / "seqproc/splitseq_bc3_map.tsv"

# Expected read structure (R2):
# NN (2bp spacer) + UMI (8bp) + BC3 (8bp) + Linker1 (30bp) + BC2 (8bp) + Linker2 (30bp) + BC1 (8bp) + ...
LINKER1 = "GTGGCCGCTGTTTCGCATCGGCGTACGACT"  # 30bp
LINKER2 = "ATCCACGTGCTTGAGAGGCCAGAGCATTCG"  # 30bp


# ============================================================================
# Helper Functions
# ============================================================================

def load_whitelist(filepath: Path) -> Dict[str, str]:
    """Load whitelist TSV file. Returns {sequence: name}."""
    whitelist = {}
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                name, seq = parts[0], parts[1]
                whitelist[seq] = name
    return whitelist


def hamming_distance(s1: str, s2: str) -> int:
    """Calculate Hamming distance between two strings."""
    if len(s1) != len(s2):
        return max(len(s1), len(s2))
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def find_closest_match(seq: str, whitelist: Dict[str, str], max_dist: int = 2) -> Tuple[str, str, int]:
    """Find closest whitelist match. Returns (name, sequence, distance)."""
    best_name, best_seq, best_dist = None, None, 100
    for wl_seq, wl_name in whitelist.items():
        if len(wl_seq) != len(seq):
            continue
        dist = hamming_distance(seq, wl_seq)
        if dist < best_dist:
            best_dist = dist
            best_name = wl_name
            best_seq = wl_seq
            if dist == 0:
                break
    return best_name, best_seq, best_dist


def find_linker_position(read: str, linker: str, max_dist: int = 3) -> Tuple[int, int]:
    """Find linker position in read with fuzzy matching. Returns (start, distance)."""
    best_pos, best_dist = -1, 100
    linker_len = len(linker)
    
    for i in range(len(read) - linker_len + 1):
        window = read[i:i + linker_len]
        dist = hamming_distance(window, linker)
        if dist < best_dist:
            best_dist = dist
            best_pos = i
            if dist == 0:
                break
    
    if best_dist <= max_dist:
        return best_pos, best_dist
    return -1, best_dist


def parse_fastq(filepath: str) -> Dict[str, Tuple[str, str]]:
    """Parse FASTQ file. Returns {read_id: (sequence, quality)}."""
    reads = {}
    with open(filepath) as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # +
            qual = f.readline().strip()
            read_id = header.strip().split()[0].replace('@', '')
            reads[read_id] = (seq, qual)
    return reads


def extract_barcodes_from_read(read: str, linker1_dist: int = 3, linker2_dist: int = 3) -> Dict:
    """
    Extract barcodes from a raw R2 read using linker anchoring.
    
    Expected structure:
    NN (2bp) + UMI (10bp) + BC3 (8bp) + Linker1 (30bp) + BC2 (8bp) + Linker2 (16bp) + BC1 (8bp) + ...
    
    Returns dict with extracted components and positions.
    """
    result = {
        'umi': None, 'bc3': None, 'bc2': None, 'bc1': None,
        'l1_pos': -1, 'l1_dist': -1, 'l2_pos': -1, 'l2_dist': -1,
        'valid': False, 'reason': ''
    }
    
    if len(read) < 82:  # Minimum length: 2 + 10 + 8 + 30 + 8 + 16 + 8 = 82
        result['reason'] = 'read_too_short'
        return result
    
    # Find Linker1
    l1_pos, l1_dist = find_linker_position(read, LINKER1, linker1_dist)
    result['l1_pos'] = l1_pos
    result['l1_dist'] = l1_dist
    
    if l1_pos < 0:
        result['reason'] = f'linker1_not_found (best_dist={l1_dist})'
        return result
    
    # Extract BC3 (8bp before Linker1)
    bc3_start = l1_pos - 8
    if bc3_start < 2:  # Need at least NN (2)
        result['reason'] = 'bc3_position_invalid'
        return result
    result['bc3'] = read[bc3_start:l1_pos]
    
    # Extract UMI (8bp before BC3)
    umi_start = bc3_start - 8
    if umi_start < 2:
        result['reason'] = 'umi_position_invalid'
        return result
    result['umi'] = read[umi_start:bc3_start]
    
    # Extract BC2 (8bp after Linker1)
    bc2_start = l1_pos + len(LINKER1)
    bc2_end = bc2_start + 8
    if bc2_end > len(read):
        result['reason'] = 'bc2_beyond_read'
        return result
    result['bc2'] = read[bc2_start:bc2_end]
    
    # Find Linker2 (should be after BC2)
    search_start = bc2_end
    l2_pos, l2_dist = find_linker_position(read[search_start:], LINKER2, linker2_dist)
    if l2_pos >= 0:
        l2_pos += search_start  # Adjust to absolute position
    result['l2_pos'] = l2_pos
    result['l2_dist'] = l2_dist
    
    if l2_pos < 0:
        result['reason'] = f'linker2_not_found (best_dist={l2_dist})'
        return result
    
    # Extract BC1 (8bp after Linker2)
    bc1_start = l2_pos + len(LINKER2)
    bc1_end = bc1_start + 8
    if bc1_end > len(read):
        result['reason'] = 'bc1_beyond_read'
        return result
    result['bc1'] = read[bc1_start:bc1_end]
    
    result['valid'] = True
    return result


# ============================================================================
# Main Analysis
# ============================================================================

def run_analysis():
    print("=" * 70)
    print("SPLiT-seq READ DISCREPANCY INVESTIGATION")
    print("=" * 70)
    
    # Load whitelists
    print("\n[1/6] Loading whitelists...")
    bc1_wl = load_whitelist(BC1_WHITELIST)
    bc2_wl = load_whitelist(BC2_WHITELIST)
    bc3_wl = load_whitelist(BC3_WHITELIST)
    print(f"  BC1: {len(bc1_wl)} sequences")
    print(f"  BC2: {len(bc2_wl)} sequences")
    print(f"  BC3: {len(bc3_wl)} sequences")
    
    # Run seqproc
    print("\n[2/6] Running seqproc (raw extraction)...")
    with tempfile.TemporaryDirectory() as tmpdir:
        sp_out1 = f"{tmpdir}/seqproc_R1.fq"
        sp_out2 = f"{tmpdir}/seqproc_R2.fq"
        
        cmd = f"{SEQPROC_BIN} --geom {CONFIG_DIR}/seqproc/splitseq_raw.geom --file1 {R1_FILE} --file2 {R2_FILE} --out1 {sp_out1} --out2 {sp_out2} --threads 4"
        start = time.time()
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=PROJECT_ROOT)
        sp_runtime = time.time() - start
        
        # Parse seqproc output
        sp_reads = parse_fastq(sp_out2) if os.path.exists(sp_out2) else {}
        print(f"  Runtime: {sp_runtime:.2f}s")
        print(f"  Output reads: {len(sp_reads):,}")
        
        # Run splitcode
        print("\n[3/6] Running splitcode...")
        sc_out1 = f"{tmpdir}/splitcode_R1.fq"
        sc_out2 = f"{tmpdir}/splitcode_R2.fq"
        sc_map = f"{tmpdir}/splitcode_mapping.txt"
        
        cmd = f"{SPLITCODE_BIN} -c {CONFIG_DIR}/splitcode/splitseq_paper.config --assign -N 2 -t 4 -m {sc_map} -o {sc_out1},{sc_out2} {R1_FILE} {R2_FILE}"
        start = time.time()
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=PROJECT_ROOT)
        sc_runtime = time.time() - start
        
        # Parse splitcode output
        sc_reads = parse_fastq(sc_out2) if os.path.exists(sc_out2) else {}
        print(f"  Runtime: {sc_runtime:.2f}s")
        print(f"  Output reads: {len(sc_reads):,}")
        
        # Parse splitcode mapping to get barcode assignments
        sc_assignments = {}
        if os.path.exists(sc_map):
            with open(sc_map) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        # Format: BARCODE_COMBO \t tag_assignments \t count
                        # We need to parse this differently
                        pass
        
        # Load raw input reads
        print("\n[4/6] Loading raw input reads...")
        raw_reads = parse_fastq(str(R2_FILE))
        print(f"  Input reads: {len(raw_reads):,}")
        
        # Analyze each raw read
        print("\n[5/6] Analyzing barcode extraction from raw reads...")
        print("  (This may take a moment...)")
        
        stats = {
            'total': len(raw_reads),
            'structure_valid': 0,
            'structure_invalid': 0,
            'bc3_exact': 0, 'bc3_dist1': 0, 'bc3_dist2': 0, 'bc3_nomatch': 0,
            'bc2_exact': 0, 'bc2_dist1': 0, 'bc2_dist2': 0, 'bc2_nomatch': 0,
            'bc1_exact': 0, 'bc1_dist1': 0, 'bc1_dist2': 0, 'bc1_nomatch': 0,
            'all_bc_valid_d1': 0,  # All 3 barcodes within dist 1
            'all_bc_valid_d2': 0,  # All 3 barcodes within dist 2
            'linker1_dist': Counter(),
            'linker2_dist': Counter(),
            'failure_reasons': Counter(),
        }
        
        # Track which reads pass each criterion
        reads_with_valid_structure = set()
        reads_with_all_bc_d1 = set()
        reads_with_all_bc_d2 = set()
        
        for i, (read_id, (seq, qual)) in enumerate(raw_reads.items()):
            if i % 100000 == 0 and i > 0:
                print(f"    Processed {i:,} reads...")
            
            # Extract barcodes using linker anchoring
            extracted = extract_barcodes_from_read(seq, linker1_dist=3, linker2_dist=3)
            
            if not extracted['valid']:
                stats['structure_invalid'] += 1
                stats['failure_reasons'][extracted['reason']] += 1
                continue
            
            stats['structure_valid'] += 1
            reads_with_valid_structure.add(read_id)
            
            # Track linker distances
            stats['linker1_dist'][extracted['l1_dist']] += 1
            stats['linker2_dist'][extracted['l2_dist']] += 1
            
            # Check BC3 against whitelist
            bc3_name, bc3_seq, bc3_dist = find_closest_match(extracted['bc3'], bc3_wl)
            if bc3_dist == 0:
                stats['bc3_exact'] += 1
            elif bc3_dist == 1:
                stats['bc3_dist1'] += 1
            elif bc3_dist == 2:
                stats['bc3_dist2'] += 1
            else:
                stats['bc3_nomatch'] += 1
            
            # Check BC2 against whitelist
            bc2_name, bc2_seq, bc2_dist = find_closest_match(extracted['bc2'], bc2_wl)
            if bc2_dist == 0:
                stats['bc2_exact'] += 1
            elif bc2_dist == 1:
                stats['bc2_dist1'] += 1
            elif bc2_dist == 2:
                stats['bc2_dist2'] += 1
            else:
                stats['bc2_nomatch'] += 1
            
            # Check BC1 against whitelist
            bc1_name, bc1_seq, bc1_dist = find_closest_match(extracted['bc1'], bc1_wl)
            if bc1_dist == 0:
                stats['bc1_exact'] += 1
            elif bc1_dist == 1:
                stats['bc1_dist1'] += 1
            elif bc1_dist == 2:
                stats['bc1_dist2'] += 1
            else:
                stats['bc1_nomatch'] += 1
            
            # Check if all barcodes are valid
            if bc3_dist <= 1 and bc2_dist <= 1 and bc1_dist <= 1:
                stats['all_bc_valid_d1'] += 1
                reads_with_all_bc_d1.add(read_id)
            
            if bc3_dist <= 2 and bc2_dist <= 2 and bc1_dist <= 2:
                stats['all_bc_valid_d2'] += 1
                reads_with_all_bc_d2.add(read_id)
        
        # Compare with tool outputs
        sp_read_ids = set(sp_reads.keys())
        sc_read_ids = set(sc_reads.keys())
        
        print("\n[6/6] Generating report...")
        
        print("\n" + "=" * 70)
        print("ANALYSIS RESULTS")
        print("=" * 70)
        
        print("\n### 1. RAW READ RECOVERY ###")
        print(f"  Total input reads:     {stats['total']:>12,}")
        print(f"  seqproc output:        {len(sp_reads):>12,} ({100*len(sp_reads)/stats['total']:.1f}%)")
        print(f"  splitcode output:      {len(sc_reads):>12,} ({100*len(sc_reads)/stats['total']:.1f}%)")
        print(f"  Difference:            {len(sc_reads) - len(sp_reads):>12,}")
        
        print("\n### 2. STRUCTURE VALIDITY (Linker Anchoring) ###")
        print(f"  Valid structure:       {stats['structure_valid']:>12,} ({100*stats['structure_valid']/stats['total']:.1f}%)")
        print(f"  Invalid structure:     {stats['structure_invalid']:>12,} ({100*stats['structure_invalid']/stats['total']:.1f}%)")
        
        print("\n  Failure reasons:")
        for reason, count in stats['failure_reasons'].most_common(10):
            print(f"    {reason}: {count:,}")
        
        print("\n  Linker1 distance distribution (of valid reads):")
        for dist in sorted(stats['linker1_dist'].keys()):
            count = stats['linker1_dist'][dist]
            print(f"    Distance {dist}: {count:,} ({100*count/stats['structure_valid']:.1f}%)")
        
        print("\n  Linker2 distance distribution (of valid reads):")
        for dist in sorted(stats['linker2_dist'].keys()):
            count = stats['linker2_dist'][dist]
            print(f"    Distance {dist}: {count:,} ({100*count/stats['structure_valid']:.1f}%)")
        
        print("\n### 3. BARCODE VALIDITY (of reads with valid structure) ###")
        valid_total = stats['structure_valid']
        
        print("\n  BC3 whitelist matching:")
        print(f"    Exact match (d=0):   {stats['bc3_exact']:>12,} ({100*stats['bc3_exact']/valid_total:.1f}%)")
        print(f"    Distance 1:          {stats['bc3_dist1']:>12,} ({100*stats['bc3_dist1']/valid_total:.1f}%)")
        print(f"    Distance 2:          {stats['bc3_dist2']:>12,} ({100*stats['bc3_dist2']/valid_total:.1f}%)")
        print(f"    No match (d>2):      {stats['bc3_nomatch']:>12,} ({100*stats['bc3_nomatch']/valid_total:.1f}%)")
        
        print("\n  BC2 whitelist matching:")
        print(f"    Exact match (d=0):   {stats['bc2_exact']:>12,} ({100*stats['bc2_exact']/valid_total:.1f}%)")
        print(f"    Distance 1:          {stats['bc2_dist1']:>12,} ({100*stats['bc2_dist1']/valid_total:.1f}%)")
        print(f"    Distance 2:          {stats['bc2_dist2']:>12,} ({100*stats['bc2_dist2']/valid_total:.1f}%)")
        print(f"    No match (d>2):      {stats['bc2_nomatch']:>12,} ({100*stats['bc2_nomatch']/valid_total:.1f}%)")
        
        print("\n  BC1 whitelist matching:")
        print(f"    Exact match (d=0):   {stats['bc1_exact']:>12,} ({100*stats['bc1_exact']/valid_total:.1f}%)")
        print(f"    Distance 1:          {stats['bc1_dist1']:>12,} ({100*stats['bc1_dist1']/valid_total:.1f}%)")
        print(f"    Distance 2:          {stats['bc1_dist2']:>12,} ({100*stats['bc1_dist2']/valid_total:.1f}%)")
        print(f"    No match (d>2):      {stats['bc1_nomatch']:>12,} ({100*stats['bc1_nomatch']/valid_total:.1f}%)")
        
        print("\n### 4. COMBINED VALIDITY ###")
        print(f"  All 3 BCs within d≤1: {stats['all_bc_valid_d1']:>12,} ({100*stats['all_bc_valid_d1']/stats['total']:.1f}%)")
        print(f"  All 3 BCs within d≤2: {stats['all_bc_valid_d2']:>12,} ({100*stats['all_bc_valid_d2']/stats['total']:.1f}%)")
        
        print("\n### 5. TOOL COMPARISON ###")
        
        # Compare seqproc output with our validity analysis
        sp_in_valid_d1 = len(sp_read_ids & reads_with_all_bc_d1)
        sp_in_valid_d2 = len(sp_read_ids & reads_with_all_bc_d2)
        sp_in_valid_struct = len(sp_read_ids & reads_with_valid_structure)
        
        print(f"\n  seqproc reads analysis:")
        print(f"    Total output:                    {len(sp_reads):>10,}")
        print(f"    With valid structure:            {sp_in_valid_struct:>10,} ({100*sp_in_valid_struct/len(sp_reads):.1f}%)")
        print(f"    With all BCs valid (d≤1):        {sp_in_valid_d1:>10,} ({100*sp_in_valid_d1/len(sp_reads):.1f}%)")
        print(f"    With all BCs valid (d≤2):        {sp_in_valid_d2:>10,} ({100*sp_in_valid_d2/len(sp_reads):.1f}%)")
        
        # Compare splitcode output
        sc_in_valid_d1 = len(sc_read_ids & reads_with_all_bc_d1)
        sc_in_valid_d2 = len(sc_read_ids & reads_with_all_bc_d2)
        sc_in_valid_struct = len(sc_read_ids & reads_with_valid_structure)
        
        print(f"\n  splitcode reads analysis:")
        print(f"    Total output:                    {len(sc_reads):>10,}")
        print(f"    With valid structure:            {sc_in_valid_struct:>10,} ({100*sc_in_valid_struct/len(sc_reads):.1f}%)")
        print(f"    With all BCs valid (d≤1):        {sc_in_valid_d1:>10,} ({100*sc_in_valid_d1/len(sc_reads):.1f}%)")
        print(f"    With all BCs valid (d≤2):        {sc_in_valid_d2:>10,} ({100*sc_in_valid_d2/len(sc_reads):.1f}%)")
        
        print("\n### 6. KEY INSIGHTS ###")
        
        # Reads only in splitcode
        sc_only = sc_read_ids - sp_read_ids
        sp_only = sp_read_ids - sc_read_ids
        both = sp_read_ids & sc_read_ids
        
        print(f"\n  Read overlap:")
        print(f"    In both tools:       {len(both):>12,}")
        print(f"    Only in seqproc:     {len(sp_only):>12,}")
        print(f"    Only in splitcode:   {len(sc_only):>12,}")
        
        # Why are splitcode-only reads not in seqproc?
        sc_only_valid_struct = len(sc_only & reads_with_valid_structure)
        sc_only_valid_d1 = len(sc_only & reads_with_all_bc_d1)
        sc_only_valid_d2 = len(sc_only & reads_with_all_bc_d2)
        
        print(f"\n  Splitcode-only reads ({len(sc_only):,}) breakdown:")
        print(f"    Have valid structure:            {sc_only_valid_struct:>10,} ({100*sc_only_valid_struct/len(sc_only):.1f}%)")
        print(f"    Have all BCs valid (d≤1):        {sc_only_valid_d1:>10,} ({100*sc_only_valid_d1/len(sc_only):.1f}%)")
        print(f"    Have all BCs valid (d≤2):        {sc_only_valid_d2:>10,} ({100*sc_only_valid_d2/len(sc_only):.1f}%)")
        
        print("\n" + "=" * 70)
        print("CONCLUSION")
        print("=" * 70)
        
        print(f"""
The read discrepancy of {len(sc_reads) - len(sp_reads):,} reads is explained by:

1. THEORETICAL MAXIMUM (structure valid): {stats['structure_valid']:,} reads
   - This is the ceiling - reads where linkers can be found

2. SEQPROC OUTPUT: {len(sp_reads):,} reads
   - Expected if using strict d≤1 filtering: {stats['all_bc_valid_d1']:,}
   - Difference from expected: {len(sp_reads) - stats['all_bc_valid_d1']:,}

3. SPLITCODE OUTPUT: {len(sc_reads):,} reads
   - Expected if using d≤2 filtering: {stats['all_bc_valid_d2']:,}
   - Difference from expected: {len(sc_reads) - stats['all_bc_valid_d2']:,}

4. EFFECTIVE YIELD (after d≤1 validity filtering):
   - seqproc high-quality reads: {sp_in_valid_d1:,}
   - splitcode high-quality reads: {sc_in_valid_d1:,}
""")


if __name__ == "__main__":
    run_analysis()
