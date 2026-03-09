#!/usr/bin/env python3
"""
Discordant Read Characterization

For each dataset, takes the unique-to-<tool> read ID files produced by
concordance_analysis.py and checks whether those reads have valid structural
features (linker positions, barcode whitelist membership).

Primary question: Are splitcode's 77K+ unique SPLiT-seq PE reads structurally
valid or false positives?
"""

import json
import os
import random
import sys
from pathlib import Path
from collections import Counter

PROJECT_ROOT = Path(__file__).parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "concordance"


# ============================================================================
# Structural validators
# ============================================================================

def hamming(s1, s2):
    if len(s1) != len(s2):
        return 99
    return sum(a != b for a, b in zip(s1, s2))


def load_whitelist(path, truncate=None):
    wl = set()
    with open(path) as f:
        for line in f:
            seq = line.strip()
            if '\t' in seq:
                parts = seq.split('\t')
                seq = parts[1] if len(parts) >= 2 else parts[0]
            if truncate and len(seq) > truncate:
                seq = seq[:truncate]
            if seq:
                wl.add(seq)
    return wl


def find_linker(read, linker, start=0, max_dist=6):
    """Find best linker position in read. Returns (position, distance)."""
    best_pos, best_dist = -1, 100
    search_end = min(len(read) - len(linker) + 1, start + 50)
    for i in range(start, max(start, search_end)):
        dist = hamming(read[i:i+len(linker)], linker)
        if dist < best_dist:
            best_dist = dist
            best_pos = i
            if dist == 0:
                break
    return best_pos, best_dist


def check_wl_match(bc, wl, max_dist=1):
    """Check if barcode matches any whitelist entry within distance."""
    if bc in wl:
        return 0
    for cand in wl:
        d = hamming(bc, cand)
        if d <= max_dist:
            return d
    return 99


# ============================================================================
# SPLiT-seq PE structural check
# ============================================================================

def analyze_splitseq_pe_reads(read_ids, r2_path, bc1_wl, bc2_wl, bc3_wl):
    """Check SPLiT-seq PE reads for structural validity.
    Structure: [NN:2][UMI:10][BC3:8][L1:30][BC2:8][L2:30][BC1:6]
    """
    L1 = "GTGGCCGCTGTTTCGCATCGGCGTACGACT"  # 30bp
    L2 = "ATCCACGTGCTTGAGAGGCCAGAGCATTCG"  # 30bp

    results = {
        "total": 0,
        "l1_found_d0": 0,
        "l1_found_d3": 0,
        "l1_found_d6": 0,
        "l1_not_found": 0,
        "l2_found_d0": 0,
        "l2_found_d3": 0,
        "l2_found_d6": 0,
        "l2_not_found": 0,
        "bc3_exact": 0,
        "bc3_d1": 0,
        "bc3_miss": 0,
        "bc2_exact": 0,
        "bc2_d1": 0,
        "bc2_miss": 0,
        "bc1_exact": 0,
        "bc1_d1": 0,
        "bc1_miss": 0,
        "fully_valid_d1": 0,
        "fully_valid_d3_linker": 0,
        "structurally_invalid": 0,
        "read_lengths": Counter(),
    }

    target_ids = set(read_ids)
    seen = 0

    with open(r2_path) as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()
            f.readline()

            rid = header.strip().split()[0].lstrip('@')
            if rid not in target_ids:
                continue

            seen += 1
            results["total"] += 1
            rlen = len(seq)
            bucket = (rlen // 10) * 10
            results["read_lengths"][bucket] += 1

            # Find L1 (expected around position 20)
            l1_pos, l1_dist = find_linker(seq, L1, start=10, max_dist=6)
            if l1_dist == 0:
                results["l1_found_d0"] += 1
            elif l1_dist <= 3:
                results["l1_found_d3"] += 1
            elif l1_dist <= 6:
                results["l1_found_d6"] += 1
            else:
                results["l1_not_found"] += 1
                results["structurally_invalid"] += 1
                continue

            # Find L2 (expected around position 58)
            l2_pos, l2_dist = find_linker(seq, L2, start=l1_pos + 30 + 8, max_dist=6)
            if l2_dist == 0:
                results["l2_found_d0"] += 1
            elif l2_dist <= 3:
                results["l2_found_d3"] += 1
            elif l2_dist <= 6:
                results["l2_found_d6"] += 1
            else:
                results["l2_not_found"] += 1
                results["structurally_invalid"] += 1
                continue

            # Extract barcodes relative to linker positions
            bc3 = seq[l1_pos-8:l1_pos] if l1_pos >= 8 else ""
            bc2 = seq[l1_pos+30:l1_pos+38]
            bc1 = seq[l2_pos+30:l2_pos+36]

            # Initialize distances (99 = no match / wrong length)
            bc3_d, bc2_d, bc1_d = 99, 99, 99

            # Check BC3
            if len(bc3) == 8:
                bc3_d = check_wl_match(bc3, bc3_wl)
                if bc3_d == 0:
                    results["bc3_exact"] += 1
                elif bc3_d <= 1:
                    results["bc3_d1"] += 1
                else:
                    results["bc3_miss"] += 1
            else:
                results["bc3_miss"] += 1

            # Check BC2
            if len(bc2) == 8:
                bc2_d = check_wl_match(bc2, bc2_wl)
                if bc2_d == 0:
                    results["bc2_exact"] += 1
                elif bc2_d <= 1:
                    results["bc2_d1"] += 1
                else:
                    results["bc2_miss"] += 1
            else:
                results["bc2_miss"] += 1

            # Check BC1
            if len(bc1) == 6:
                bc1_d = check_wl_match(bc1, bc1_wl)
                if bc1_d == 0:
                    results["bc1_exact"] += 1
                elif bc1_d <= 1:
                    results["bc1_d1"] += 1
                else:
                    results["bc1_miss"] += 1
            else:
                results["bc1_miss"] += 1

            # Fully valid: both linkers within d<=3 AND all BCs within d<=1
            linker_ok = l1_dist <= 3 and l2_dist <= 3
            bc_ok = (bc3_d <= 1 and bc2_d <= 1 and bc1_d <= 1)

            if linker_ok and bc_ok:
                results["fully_valid_d1"] += 1
            elif l1_dist <= 6 and l2_dist <= 6 and bc_ok:
                results["fully_valid_d3_linker"] += 1

            if seen >= len(target_ids):
                break

    return results


def print_analysis(label, results, total_reads):
    """Pretty-print structural analysis results."""
    t = results["total"]
    if t == 0:
        print(f"  {label}: No reads to analyze")
        return

    print(f"\n  {label} ({t:,} reads)")
    print(f"  {'─'*50}")

    # Linker 1
    l1_any = results["l1_found_d0"] + results["l1_found_d3"] + results["l1_found_d6"]
    print(f"  Linker1: d=0: {results['l1_found_d0']:,} ({results['l1_found_d0']/t*100:.1f}%) | "
          f"d<=3: {results['l1_found_d3']:,} | d<=6: {results['l1_found_d6']:,} | "
          f"not found: {results['l1_not_found']:,}")

    # Linker 2
    l2_any = results["l2_found_d0"] + results["l2_found_d3"] + results["l2_found_d6"]
    print(f"  Linker2: d=0: {results['l2_found_d0']:,} ({results['l2_found_d0']/t*100:.1f}%) | "
          f"d<=3: {results['l2_found_d3']:,} | d<=6: {results['l2_found_d6']:,} | "
          f"not found: {results['l2_not_found']:,}")

    # Barcodes
    print(f"  BC3: exact={results['bc3_exact']:,} d1={results['bc3_d1']:,} miss={results['bc3_miss']:,}")
    print(f"  BC2: exact={results['bc2_exact']:,} d1={results['bc2_d1']:,} miss={results['bc2_miss']:,}")
    print(f"  BC1: exact={results['bc1_exact']:,} d1={results['bc1_d1']:,} miss={results['bc1_miss']:,}")

    # Summary
    print(f"  Fully valid (linker d<=3, BC d<=1): {results['fully_valid_d1']:,} ({results['fully_valid_d1']/t*100:.1f}%)")
    print(f"  Valid with relaxed linker (d<=6):    {results['fully_valid_d3_linker']:,}")
    print(f"  Structurally invalid (no linkers):   {results['structurally_invalid']:,} ({results['structurally_invalid']/t*100:.1f}%)")


def main():
    print("=" * 70)
    print("DISCORDANT READ CHARACTERIZATION")
    print("=" * 70)

    # Load SPLiT-seq PE whitelists
    bc23_wl = load_whitelist(PROJECT_ROOT / "configs/seqproc/splitseq_bc23_whitelist.txt")
    bc1_wl = load_whitelist(PROJECT_ROOT / "configs/seqproc/splitseq_bc1_whitelist_6bp.txt")
    print(f"Loaded whitelists: BC2/3={len(bc23_wl)} entries, BC1={len(bc1_wl)} entries (6bp)")

    r2_path = PROJECT_ROOT / "data/SRR6750041_1M_R2.fastq"
    if not r2_path.exists():
        print(f"[ERROR] R2 data not found: {r2_path}")
        sys.exit(1)

    # Load discordant read ID sets
    ds_dir = RESULTS_DIR / "splitseq_pe"
    if not ds_dir.exists():
        print(f"[ERROR] Run concordance_analysis.py --datasets splitseq_pe first")
        sys.exit(1)

    analyses = {}

    # 1. Splitcode-unique reads (the big question: are these valid?)
    sc_unique_file = ds_dir / "unique_to_splitcode.txt"
    if sc_unique_file.exists():
        with open(sc_unique_file) as f:
            sc_unique_ids = {line.strip() for line in f}
        print(f"\nAnalyzing {len(sc_unique_ids):,} splitcode-unique reads...")
        analyses["splitcode_unique"] = analyze_splitseq_pe_reads(
            sc_unique_ids, str(r2_path), bc1_wl, bc23_wl, bc23_wl)
    else:
        print(f"[SKIP] No splitcode-unique IDs found")

    # 2. Seqproc-unique reads (reads seqproc finds that no other tool does)
    sp_unique_file = ds_dir / "unique_to_seqproc.txt"
    if sp_unique_file.exists():
        with open(sp_unique_file) as f:
            sp_unique_ids = {line.strip() for line in f}
        print(f"\nAnalyzing {len(sp_unique_ids):,} seqproc-unique reads...")
        analyses["seqproc_unique"] = analyze_splitseq_pe_reads(
            sp_unique_ids, str(r2_path), bc1_wl, bc23_wl, bc23_wl)
    else:
        print(f"[SKIP] No seqproc-unique IDs found")

    # 3. Matchbox-unique reads
    mb_unique_file = ds_dir / "unique_to_matchbox.txt"
    if mb_unique_file.exists():
        with open(mb_unique_file) as f:
            mb_unique_ids = {line.strip() for line in f}
        print(f"\nAnalyzing {len(mb_unique_ids):,} matchbox-unique reads...")
        if mb_unique_ids:
            analyses["matchbox_unique"] = analyze_splitseq_pe_reads(
                mb_unique_ids, str(r2_path), bc1_wl, bc23_wl, bc23_wl)
    else:
        print(f"[SKIP] No matchbox-unique IDs found")

    # 4. Also check splitcode reads that are NOT in seqproc
    # (broader set than "unique to splitcode" -- includes those shared with matchbox)
    sc_ids_file = ds_dir / "splitcode_ids.txt"
    sp_ids_file = ds_dir / "seqproc_edit_ids.txt"
    mb_ids_file = ds_dir / "matchbox_ids.txt"

    # Assert required ID files exist before cross-set analysis
    for required_file in [sc_ids_file, sp_ids_file, mb_ids_file]:
        if not required_file.exists():
            print(f"[ERROR] Required ID file missing: {required_file}")
            print(f"Run concordance_analysis.py --datasets splitseq_pe first.")
            sys.exit(1)

    with open(sc_ids_file) as f:
        sc_ids = {line.strip() for line in f}
    with open(sp_ids_file) as f:
        sp_ids = {line.strip() for line in f}
    with open(mb_ids_file) as f:
        mb_ids = {line.strip() for line in f}

    sc_not_sp = sc_ids - sp_ids
    if sc_not_sp:
        print(f"\nAnalyzing {len(sc_not_sp):,} splitcode-but-not-seqproc reads...")
        analyses["splitcode_not_seqproc"] = analyze_splitseq_pe_reads(
            sc_not_sp, str(r2_path), bc1_wl, bc23_wl, bc23_wl)

    # 5. Also sample some consensus reads as a control
    consensus = sp_ids & mb_ids & sc_ids
    # Sample 10K from consensus for control
    random.seed(42)
    consensus_sample = set(random.sample(sorted(consensus), min(10000, len(consensus))))
    print(f"\nAnalyzing {len(consensus_sample):,} consensus reads (control sample)...")
    analyses["consensus_control"] = analyze_splitseq_pe_reads(
        consensus_sample, str(r2_path), bc1_wl, bc23_wl, bc23_wl)

    # Print results
    print(f"\n{'='*70}")
    print("RESULTS")
    print(f"{'='*70}")

    for label, res in analyses.items():
        print_analysis(label, res, 1_000_000)

    # Save results as JSON
    # Convert Counter objects for JSON serialization
    json_analyses = {}
    for label, res in analyses.items():
        r = dict(res)
        r["read_lengths"] = dict(r["read_lengths"])
        json_analyses[label] = r

    out_file = ds_dir / "discordant_analysis.json"
    with open(out_file, 'w') as f:
        json.dump(json_analyses, f, indent=2)
    print(f"\nSaved: {out_file}")

    # Key conclusion
    if "splitcode_unique" in analyses:
        sc = analyses["splitcode_unique"]
        t = sc["total"]
        valid = sc["fully_valid_d1"]
        invalid = sc["structurally_invalid"]
        print(f"\n{'='*70}")
        print("KEY FINDING: splitcode-unique SPLiT-seq PE reads")
        print(f"{'='*70}")
        print(f"  Total unique to splitcode: {t:,}")
        print(f"  Structurally valid (linker d<=3, BC d<=1): {valid:,} ({valid/max(t,1)*100:.1f}%)")
        print(f"  Structurally invalid (no linkers): {invalid:,} ({invalid/max(t,1)*100:.1f}%)")
        if valid / max(t, 1) < 0.5:
            print(f"  CONCLUSION: Majority of splitcode-unique reads are FALSE POSITIVES")
        else:
            print(f"  CONCLUSION: Majority of splitcode-unique reads appear structurally valid")


if __name__ == "__main__":
    main()
