#!/usr/bin/env python3
"""
Quick Jaccard Index Analysis for Sci-Seq 3 Dataset
Compares read ID sets between seqproc, matchbox, and splitcode outputs.
"""

import subprocess
import tempfile
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"

# Tool binaries
SEQPROC_BIN = os.environ.get("SEQPROC_BIN", str(PROJECT_ROOT.parent / "seqproc/target/release/seqproc"))
MATCHBOX_BIN = os.environ.get("MATCHBOX_BIN", str(PROJECT_ROOT.parent / "matchbox/target/release/matchbox"))
SPLITCODE_BIN = os.environ.get("SPLITCODE_BIN", str(PROJECT_ROOT.parent / "splitcode/build/src/splitcode"))

# Sciseq dataset
R1 = DATA_DIR / "SRR7827254_1M_1.fastq"
R2 = DATA_DIR / "SRR7827254_1M_2.fastq"

# Configs
SEQPROC_GEOM = PROJECT_ROOT / "configs/seqproc/sciseq3.geom"
MATCHBOX_CONFIG = PROJECT_ROOT / "configs/matchbox/sciseq3.mb"
SPLITCODE_CONFIG = PROJECT_ROOT / "configs/splitcode/sciseq3.config"


def extract_read_ids(fastq_path):
    """Extract read IDs from FASTQ file."""
    ids = set()
    with open(fastq_path, 'r') as f:
        while True:
            header = f.readline()
            if not header:
                break
            # Extract ID from header (everything before first space, without @)
            read_id = header.strip().split()[0].replace('@', '')
            ids.add(read_id)
            f.readline()  # seq
            f.readline()  # +
            f.readline()  # qual
    return ids


def jaccard(set_a, set_b):
    """Compute Jaccard index."""
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    return intersection / union if union > 0 else 0.0


def run_seqproc(tmpdir):
    """Run seqproc and return output read IDs."""
    out1 = f"{tmpdir}/seqproc_R1.fq"
    out2 = f"{tmpdir}/seqproc_R2.fq"
    cmd = f"{SEQPROC_BIN} --geom {SEQPROC_GEOM} --file1 {R1} --file2 {R2} --out1 {out1} --out2 {out2} --threads 4"
    print(f"  Running seqproc...")
    subprocess.run(cmd, shell=True, capture_output=True, cwd=PROJECT_ROOT)
    return extract_read_ids(out2)


def run_matchbox(tmpdir):
    """Run matchbox and return output read IDs."""
    out_tsv = f"{tmpdir}/matchbox_out.tsv"
    cmd = f'{MATCHBOX_BIN} -e 0.2 -t 4 -s {MATCHBOX_CONFIG} {R1} -p {R2} > {out_tsv}'
    print(f"  Running matchbox...")
    subprocess.run(cmd, shell=True, capture_output=True, cwd=PROJECT_ROOT)
    
    # Matchbox generates mb_r1.fq, mb_r2.fq in PROJECT_ROOT
    mb_r2 = PROJECT_ROOT / "mb_r2.fq"
    if mb_r2.exists():
        ids = extract_read_ids(mb_r2)
        # Clean up
        (PROJECT_ROOT / "mb_r1.fq").unlink(missing_ok=True)
        mb_r2.unlink(missing_ok=True)
        return ids
    return set()


def run_splitcode(tmpdir):
    """Run splitcode and return output read IDs."""
    out1 = f"{tmpdir}/splitcode_R1.fq"
    out2 = f"{tmpdir}/splitcode_R2.fq"
    mapping = f"{tmpdir}/splitcode_mapping.txt"
    cmd = f"{SPLITCODE_BIN} -c {SPLITCODE_CONFIG} --assign -N 2 -t 4 -m {mapping} -o {out1},{out2} {R1} {R2}"
    print(f"  Running splitcode...")
    subprocess.run(cmd, shell=True, capture_output=True, cwd=PROJECT_ROOT)
    return extract_read_ids(out2)


def main():
    print("=" * 60)
    print("Sci-Seq 3 Jaccard Index Analysis")
    print("=" * 60)
    print(f"Input R1: {R1}")
    print(f"Input R2: {R2}")
    print()
    
    # Check inputs exist
    if not R1.exists() or not R2.exists():
        print(f"ERROR: Input files not found!")
        return
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Run all tools
        seqproc_ids = run_seqproc(tmpdir)
        matchbox_ids = run_matchbox(tmpdir)
        splitcode_ids = run_splitcode(tmpdir)
        
        print()
        print("=" * 60)
        print("RESULTS")
        print("=" * 60)
        
        print(f"\n## Read Counts")
        print(f"  seqproc:   {len(seqproc_ids):,} reads")
        print(f"  matchbox:  {len(matchbox_ids):,} reads")
        print(f"  splitcode: {len(splitcode_ids):,} reads")
        
        print(f"\n## Pairwise Intersections")
        sp_mb = seqproc_ids & matchbox_ids
        sp_sc = seqproc_ids & splitcode_ids
        mb_sc = matchbox_ids & splitcode_ids
        all_three = seqproc_ids & matchbox_ids & splitcode_ids
        
        print(f"  seqproc ∩ matchbox:  {len(sp_mb):,}")
        print(f"  seqproc ∩ splitcode: {len(sp_sc):,}")
        print(f"  matchbox ∩ splitcode: {len(mb_sc):,}")
        print(f"  All three:           {len(all_three):,}")
        
        print(f"\n## Jaccard Indices")
        j_sp_mb = jaccard(seqproc_ids, matchbox_ids)
        j_sp_sc = jaccard(seqproc_ids, splitcode_ids)
        j_mb_sc = jaccard(matchbox_ids, splitcode_ids)
        
        print(f"  seqproc vs matchbox:  {j_sp_mb:.4f}")
        print(f"  seqproc vs splitcode: {j_sp_sc:.4f}")
        print(f"  matchbox vs splitcode: {j_mb_sc:.4f}")
        
        print(f"\n## Unique to Each Tool")
        only_seqproc = seqproc_ids - matchbox_ids - splitcode_ids
        only_matchbox = matchbox_ids - seqproc_ids - splitcode_ids
        only_splitcode = splitcode_ids - seqproc_ids - matchbox_ids
        
        print(f"  Only seqproc:   {len(only_seqproc):,}")
        print(f"  Only matchbox:  {len(only_matchbox):,}")
        print(f"  Only splitcode: {len(only_splitcode):,}")


if __name__ == "__main__":
    main()
