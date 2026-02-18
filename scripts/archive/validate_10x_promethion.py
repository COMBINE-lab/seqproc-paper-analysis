#!/usr/bin/env python3
"""
10x PromethION Barcode Validity Analysis

This script compares seqproc and matchbox on 10x PromethION long-read data,
validating extracted barcodes against the official 10x Genomics whitelist.

VALIDATION METHOD:
1. Run seqproc with forward primer geometry (finds reads with primer in forward orientation)
2. Run seqproc with reverse primer geometry (finds reads with primer in reverse complement)
3. Merge seqproc results, keeping unique read IDs
4. Run matchbox on the same input
5. Extract 16bp barcodes from each tool's output
6. Validate barcodes against 10x Genomics 3M whitelist (exact match)
7. Compute overlap statistics and validity rates for each subset

The 10x whitelist contains ~6.8M barcodes (including reverse complements).
A barcode is considered "valid" if it exactly matches one of these sequences.
"""

import os
import sys
import gzip
import tempfile
import subprocess
from pathlib import Path
from collections import defaultdict

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "10x"
INPUT_FASTQ = DATA_DIR / "ERR9958135_1M.fastq"
WHITELIST_PATH = Path("/home/ubuntu/3M-february-2018.txt.gz")

# Tool paths
SEQPROC_BIN = PROJECT_ROOT.parent / "seqproc" / "target" / "release" / "seqproc"
MATCHBOX_BIN = PROJECT_ROOT.parent / "matchbox" / "target" / "release" / "matchbox"
MATCHBOX_CONFIG = PROJECT_ROOT / "configs" / "matchbox" / "10x_longread.mb"

# Primer sequences
PRIMER_FWD = "CTACACGACGCTCTTCCGATCT"  # 22bp
PRIMER_REV = "AGATCGGAAGAGCGTCGTGTAG"  # reverse complement


def load_whitelist(path: Path) -> set:
    """Load 10x barcode whitelist from gzipped file."""
    whitelist = set()
    with gzip.open(path, 'rt') as f:
        for line in f:
            whitelist.add(line.strip())
    return whitelist


def create_seqproc_geometry(outdir: Path, direction: str) -> Path:
    """
    Create seqproc geometry file for 10x long-read extraction.
    
    Forward: primer -> barcode -> UMI -> rest
    Reverse: rest -> UMI -> barcode -> primer (reverse complement)
    """
    geom_path = outdir / f"10x_lr_{direction}.geom"
    
    if direction == "fwd":
        content = f"""primer = anchor_relative(hamming(f[{PRIMER_FWD}], 3))
bc = b[16]
umi = u[12]
rest = r:
1{{<primer><bc><umi><rest>}}
-> 1{{<bc><umi>}}
"""
    else:  # reverse
        content = f"""primer = anchor_relative(hamming(f[{PRIMER_REV}], 3))
bc = revcomp(b[16])
umi = revcomp(u[12])
rest = r:
1{{<rest><umi><bc><primer>}}
-> 1{{<bc><umi>}}
"""
    
    with open(geom_path, 'w') as f:
        f.write(content)
    
    return geom_path


def run_seqproc(input_fastq: Path, geom_path: Path, output_fastq: Path) -> int:
    """Run seqproc and return the number of output reads."""
    cmd = [
        str(SEQPROC_BIN),
        "--geom", str(geom_path),
        "--file1", str(input_fastq),
        "--out1", str(output_fastq),
        "--threads", "4"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Seqproc error: {result.stderr}", file=sys.stderr)
        return 0
    
    # Count reads
    count = 0
    with open(output_fastq) as f:
        for line in f:
            if line.startswith('@'):
                count += 1
    return count


def run_matchbox(input_fastq: Path, config_path: Path, output_tsv: Path) -> int:
    """Run matchbox and return the number of output reads."""
    cmd = [
        str(MATCHBOX_BIN),
        "-e", "0.2",
        "-t", "4",
        "-s", str(config_path),
        str(input_fastq)
    ]
    
    with open(output_tsv, 'w') as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
    
    if result.returncode != 0:
        print(f"Matchbox error: {result.stderr}", file=sys.stderr)
        return 0
    
    # Count lines
    count = 0
    with open(output_tsv) as f:
        for _ in f:
            count += 1
    return count


def extract_seqproc_barcodes(fastq_paths: list) -> dict:
    """
    Extract barcodes from seqproc output FASTQ files.
    
    Seqproc outputs: BC (16bp) + UMI (12bp) = 28bp sequence
    We extract the first 16bp as the barcode.
    
    Returns: {read_id: barcode}
    """
    barcodes = {}
    
    for fastq_path in fastq_paths:
        with open(fastq_path) as f:
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip()
                f.readline()  # +
                f.readline()  # quality
                
                # Parse read ID (remove @ and any trailing info)
                read_id = header.split()[0][1:]
                
                # Extract barcode (first 16bp)
                bc = seq[:16]
                
                # Keep first occurrence if duplicate read IDs
                if read_id not in barcodes:
                    barcodes[read_id] = bc
    
    return barcodes


def extract_matchbox_barcodes(tsv_path: Path) -> dict:
    """
    Extract barcodes from matchbox TSV output.
    
    Matchbox outputs tab-separated: read_id, barcode, [umi]
    Barcode is 16bp.
    
    Returns: {read_id: barcode}
    """
    barcodes = {}
    
    with open(tsv_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                read_id = parts[0]
                bc = parts[1][:16] if len(parts[1]) >= 16 else parts[1]
                
                if read_id not in barcodes:
                    barcodes[read_id] = bc
    
    return barcodes


def validate_barcodes(barcodes: dict, whitelist: set) -> tuple:
    """
    Validate barcodes against whitelist.
    
    Returns: (total_count, valid_count)
    """
    total = len(barcodes)
    valid = sum(1 for bc in barcodes.values() if bc in whitelist)
    return total, valid


def main():
    print("=" * 70)
    print("10x PromethION Barcode Validity Analysis")
    print("=" * 70)
    
    # Check prerequisites
    if not INPUT_FASTQ.exists():
        print(f"ERROR: Input file not found: {INPUT_FASTQ}")
        sys.exit(1)
    
    if not WHITELIST_PATH.exists():
        print(f"ERROR: Whitelist not found: {WHITELIST_PATH}")
        sys.exit(1)
    
    if not SEQPROC_BIN.exists():
        print(f"ERROR: Seqproc binary not found: {SEQPROC_BIN}")
        sys.exit(1)
    
    if not MATCHBOX_BIN.exists():
        print(f"ERROR: Matchbox binary not found: {MATCHBOX_BIN}")
        sys.exit(1)
    
    # Load whitelist
    print(f"\n[1] Loading 10x whitelist from {WHITELIST_PATH}...")
    whitelist = load_whitelist(WHITELIST_PATH)
    print(f"    Loaded {len(whitelist):,} barcodes")
    
    # Create temp directory
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Run seqproc (forward)
        print(f"\n[2] Running seqproc (forward orientation)...")
        geom_fwd = create_seqproc_geometry(tmpdir, "fwd")
        seqproc_fwd_out = tmpdir / "seqproc_fwd.fq"
        fwd_count = run_seqproc(INPUT_FASTQ, geom_fwd, seqproc_fwd_out)
        print(f"    Found {fwd_count:,} reads with forward primer")
        
        # Run seqproc (reverse)
        print(f"\n[3] Running seqproc (reverse orientation)...")
        geom_rev = create_seqproc_geometry(tmpdir, "rev")
        seqproc_rev_out = tmpdir / "seqproc_rev.fq"
        rev_count = run_seqproc(INPUT_FASTQ, geom_rev, seqproc_rev_out)
        print(f"    Found {rev_count:,} reads with reverse primer")
        
        # Run matchbox
        print(f"\n[4] Running matchbox...")
        matchbox_out = tmpdir / "matchbox.tsv"
        mb_count = run_matchbox(INPUT_FASTQ, MATCHBOX_CONFIG, matchbox_out)
        print(f"    Found {mb_count:,} reads")
        
        # Extract barcodes
        print(f"\n[5] Extracting barcodes...")
        seqproc_barcodes = extract_seqproc_barcodes([seqproc_fwd_out, seqproc_rev_out])
        matchbox_barcodes = extract_matchbox_barcodes(matchbox_out)
        print(f"    Seqproc unique reads: {len(seqproc_barcodes):,}")
        print(f"    Matchbox unique reads: {len(matchbox_barcodes):,}")
        
        # Validate barcodes
        print(f"\n[6] Validating barcodes against whitelist...")
        seq_total, seq_valid = validate_barcodes(seqproc_barcodes, whitelist)
        mb_total, mb_valid = validate_barcodes(matchbox_barcodes, whitelist)
        
        # Compute overlap
        seqproc_ids = set(seqproc_barcodes.keys())
        matchbox_ids = set(matchbox_barcodes.keys())
        
        intersection = seqproc_ids & matchbox_ids
        seqproc_only = seqproc_ids - matchbox_ids
        matchbox_only = matchbox_ids - seqproc_ids
        
        # Validate subsets
        int_seq_valid = sum(1 for rid in intersection if seqproc_barcodes[rid] in whitelist)
        int_mb_valid = sum(1 for rid in intersection if matchbox_barcodes[rid] in whitelist)
        seq_only_valid = sum(1 for rid in seqproc_only if seqproc_barcodes[rid] in whitelist)
        mb_only_valid = sum(1 for rid in matchbox_only if matchbox_barcodes[rid] in whitelist)
        
        # Barcode agreement in intersection
        agree_count = sum(1 for rid in intersection 
                        if seqproc_barcodes[rid] == matchbox_barcodes[rid])
        
        # Print results
        print("\n" + "=" * 70)
        print("RESULTS")
        print("=" * 70)
        
        print(f"\n### OVERALL RECOVERY ###")
        print(f"Seqproc:  {seq_total:>10,} reads ({100*seq_total/1_000_000:.1f}% of input)")
        print(f"Matchbox: {mb_total:>10,} reads ({100*mb_total/1_000_000:.1f}% of input)")
        
        print(f"\n### OVERALL VALIDITY (exact whitelist match) ###")
        print(f"Seqproc:  {seq_valid:>10,} valid ({100*seq_valid/seq_total:.1f}%)")
        print(f"Matchbox: {mb_valid:>10,} valid ({100*mb_valid/mb_total:.1f}%)")
        
        print(f"\n### OVERLAP ANALYSIS ###")
        print(f"Intersection:   {len(intersection):>10,} reads")
        print(f"Seqproc-only:   {len(seqproc_only):>10,} reads")
        print(f"Matchbox-only:  {len(matchbox_only):>10,} reads")
        
        print(f"\n### VALIDITY BY SUBSET ###")
        if len(intersection) > 0:
            print(f"Intersection (seqproc BC):  {int_seq_valid:>10,} valid ({100*int_seq_valid/len(intersection):.1f}%)")
            print(f"Intersection (matchbox BC): {int_mb_valid:>10,} valid ({100*int_mb_valid/len(intersection):.1f}%)")
        if len(seqproc_only) > 0:
            print(f"Seqproc-only:               {seq_only_valid:>10,} valid ({100*seq_only_valid/len(seqproc_only):.1f}%)")
        if len(matchbox_only) > 0:
            print(f"Matchbox-only:              {mb_only_valid:>10,} valid ({100*mb_only_valid/len(matchbox_only):.1f}%)")
        
        print(f"\n### BARCODE AGREEMENT ###")
        if len(intersection) > 0:
            print(f"Reads where both tools extract same barcode: {agree_count:,} ({100*agree_count/len(intersection):.1f}%)")
        
        print("\n" + "=" * 70)
        print("INTERPRETATION")
        print("=" * 70)
        print(f"""
- Seqproc recovers {seq_total:,} reads with {100*seq_valid/seq_total:.1f}% validity
- Matchbox recovers {mb_total:,} reads with {100*mb_valid/mb_total:.1f}% validity
- Matchbox finds {len(matchbox_only):,} additional reads not found by seqproc
- Of those additional reads, {mb_only_valid:,} ({100*mb_only_valid/len(matchbox_only):.1f}%) have valid barcodes
- When both tools find the same read, they agree on barcode {100*agree_count/len(intersection):.1f}% of the time

This represents a PRECISION vs RECALL tradeoff:
- Seqproc: Higher precision ({100*seq_valid/seq_total:.1f}%), lower recall ({100*seq_total/1_000_000:.1f}%)
- Matchbox: Lower precision ({100*mb_valid/mb_total:.1f}%), higher recall ({100*mb_total/1_000_000:.1f}%)
""")


if __name__ == "__main__":
    main()
