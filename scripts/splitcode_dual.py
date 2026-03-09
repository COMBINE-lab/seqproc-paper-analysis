#!/usr/bin/env python3
"""
Splitcode Dual-Orientation Analysis for LR-SPLiT-seq

Splitcode does not natively support dual-orientation matching. This script
simulates it by:
  1. Running splitcode on the original (forward) reads (uses cached IDs)
  2. Reverse-complementing all reads
  3. Running splitcode on the RC'd reads
  4. Combining (union) the read ID sets from both runs
  5. Reporting combined recovery and concordance vs seqproc/matchbox

This is a SUPPLEMENTARY analysis to show what splitcode could achieve if it
had orientation support. The runtime is effectively doubled.

Usage:
    python3 scripts/splitcode_dual.py --threads 4
    python3 scripts/splitcode_dual.py --threads 4 --skip-rc-run
"""

import json
import os
import subprocess
import time
import argparse
from pathlib import Path
from typing import Set, Tuple

PROJECT_ROOT = Path(__file__).parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "concordance" / "lr_splitseq"

SPLITCODE_BIN = os.environ.get(
    "SPLITCODE_BIN",
    str(PROJECT_ROOT.parent / "splitcode/build/src/splitcode")
)

LR_FASTQ = PROJECT_ROOT / "data" / "SRR13948564_1M.fastq"
SPLITCODE_CONFIG = PROJECT_ROOT / "configs" / "splitcode" / "splitseq_singleend.config"

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]


def rc_fastq(input_path: Path, output_path: Path) -> int:
    """Reverse-complement all reads in a FASTQ file. Returns read count."""
    count = 0
    with open(input_path) as fin, open(output_path, 'w') as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline().strip()
            plus = fin.readline()
            qual = fin.readline().strip()

            fout.write(header)
            fout.write(reverse_complement(seq) + '\n')
            fout.write(plus)
            fout.write(qual[::-1] + '\n')
            count += 1
    return count


def extract_fastq_ids(filepath: str) -> Set[str]:
    """Extract read IDs from a FASTQ file."""
    ids = set()
    if not os.path.exists(filepath):
        return ids
    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                rid = line.strip().split()[0].lstrip('@')
                ids.add(rid)
    return ids


def load_id_file(filepath: Path) -> Set[str]:
    """Load read IDs from a text file (one per line)."""
    if not filepath.exists():
        return set()
    with open(filepath) as f:
        return {line.strip() for line in f if line.strip()}


def jaccard(a: Set[str], b: Set[str]) -> float:
    """Compute Jaccard index."""
    if not a and not b:
        return 1.0
    union = len(a | b)
    return len(a & b) / union if union else 1.0


def run_cmd(cmd: str, cwd=None) -> Tuple[float, float, int, str]:
    """Run command, return (runtime_s, peak_mem_mb, exit_code, stderr)."""
    time_cmd = f"/usr/bin/time -v {cmd}"
    start = time.time()
    result = subprocess.run(time_cmd, shell=True, capture_output=True, text=True, cwd=cwd)
    runtime = time.time() - start

    peak_mem_kb = 0
    tool_exit_code = result.returncode
    for line in result.stderr.split('\n'):
        if 'Maximum resident set size' in line:
            peak_mem_kb = int(line.split(':')[1].strip())
        elif 'Exit status' in line:
            tool_exit_code = int(line.split(':')[1].strip())

    return runtime, peak_mem_kb / 1024, tool_exit_code, result.stderr


def main():
    parser = argparse.ArgumentParser(
        description='Splitcode dual-orientation analysis for LR-SPLiT-seq')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--skip-rc-run', action='store_true',
                        help='Skip RC splitcode run, use cached IDs')
    args = parser.parse_args()

    print("=" * 70)
    print("SUPPLEMENTARY: SPLITCODE DUAL-ORIENTATION (LR-SPLiT-seq)")
    print("=" * 70)
    print("NOTE: This is a simulated dual-orientation run. Splitcode does not")
    print("natively support orientation-aware matching. Runtime is 2x.")
    print()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    if not LR_FASTQ.exists():
        print(f"[ERROR] Input FASTQ not found: {LR_FASTQ}")
        return

    # Step 1: Load forward splitcode IDs (from concordance_analysis.py cache)
    fw_ids_file = RESULTS_DIR / "splitcode_ids.txt"
    if not fw_ids_file.exists():
        print(f"[ERROR] Forward splitcode IDs not found: {fw_ids_file}")
        print("Run: python3 scripts/concordance_analysis.py --datasets lr_splitseq")
        return

    fw_ids = load_id_file(fw_ids_file)
    print(f"Forward splitcode: {len(fw_ids):,} reads (cached)")

    # Step 2: Reverse complement the FASTQ
    rc_fastq_path = RESULTS_DIR / "SRR13948564_1M_RC.fastq"
    if not rc_fastq_path.exists():
        print(f"Reverse complementing {LR_FASTQ.name}...", end=" ", flush=True)
        count = rc_fastq(LR_FASTQ, rc_fastq_path)
        print(f"{count:,} reads")
    else:
        print(f"RC FASTQ exists: {rc_fastq_path.name}")

    # Step 3: Run splitcode on RC reads
    rc_ids_file = RESULTS_DIR / "splitcode_rc_ids.txt"
    rc_runtime, rc_mem = 0, 0

    if args.skip_rc_run and rc_ids_file.exists():
        print(f"RC splitcode: loading cached IDs...")
        rc_ids = load_id_file(rc_ids_file)
    else:
        rc_out = RESULTS_DIR / "splitcode_rc_out.fq"
        rc_mapping = RESULTS_DIR / "splitcode_rc_mapping.txt"

        cmd = (f"{SPLITCODE_BIN} -c {SPLITCODE_CONFIG} "
               f"--assign -N 1 -t {args.threads} -m {rc_mapping} "
               f"-o {rc_out} {rc_fastq_path}")

        print(f"Running splitcode on RC reads...", end=" ", flush=True)
        rc_runtime, rc_mem, rc_code, stderr = run_cmd(cmd, PROJECT_ROOT)

        if rc_code != 0:
            print(f"\n[ERROR] splitcode RC failed (rc={rc_code})")
            print(stderr[-500:])
            return

        rc_ids = extract_fastq_ids(str(rc_out))

        # Save RC IDs
        with open(rc_ids_file, 'w') as f:
            for rid in sorted(rc_ids):
                f.write(rid + '\n')

        print(f"{len(rc_ids):,} reads in {rc_runtime:.1f}s, {rc_mem:.0f}MB")

    # Step 4: Combine forward + RC IDs
    combined_ids = fw_ids | rc_ids
    both_ori = fw_ids & rc_ids
    fw_only = fw_ids - rc_ids
    rc_only = rc_ids - fw_ids

    print(f"\n{'='*60}")
    print("SPLITCODE DUAL-ORIENTATION RESULTS")
    print(f"{'='*60}")
    print(f"  Forward only:     {len(fw_ids):>10,} reads")
    print(f"  RC only:          {len(rc_ids):>10,} reads")
    print(f"  Both orientations:{len(both_ori):>10,} reads")
    print(f"  Forward-unique:   {len(fw_only):>10,} reads")
    print(f"  RC-unique:        {len(rc_only):>10,} reads")
    print(f"  Combined (union): {len(combined_ids):>10,} reads")
    print(f"  Combined recovery:{len(combined_ids)/1_000_000*100:>9.1f}%")

    # Step 5: Concordance with seqproc and matchbox
    print(f"\n--- Concordance with other tools ---")

    sp_ids = load_id_file(RESULTS_DIR / "seqproc_edit_ids.txt")
    mb_ids = load_id_file(RESULTS_DIR / "matchbox_ids.txt")

    if sp_ids:
        j = jaccard(combined_ids, sp_ids)
        both = len(combined_ids & sp_ids)
        sc_only = len(combined_ids - sp_ids)
        sp_only = len(sp_ids - combined_ids)
        print(f"  splitcode_dual vs seqproc:  J={j:.4f}  "
              f"(both={both:,}, sc_dual_only={sc_only:,}, sp_only={sp_only:,})")

    if mb_ids:
        j = jaccard(combined_ids, mb_ids)
        both = len(combined_ids & mb_ids)
        sc_only = len(combined_ids - mb_ids)
        mb_only = len(mb_ids - combined_ids)
        print(f"  splitcode_dual vs matchbox: J={j:.4f}  "
              f"(both={both:,}, sc_dual_only={sc_only:,}, mb_only={mb_only:,})")

    # Compare improvement over forward-only
    print(f"\n--- Improvement over forward-only ---")
    gain = len(combined_ids) - len(fw_ids)
    print(f"  Forward:  {len(fw_ids):>10,} ({len(fw_ids)/1_000_000*100:.1f}%)")
    print(f"  Combined: {len(combined_ids):>10,} ({len(combined_ids)/1_000_000*100:.1f}%)")
    print(f"  Gain:     +{gain:>9,} (+{gain/max(len(fw_ids),1)*100:.1f}%)")

    # Still compare to seqproc ann+edit
    if sp_ids:
        print(f"\n--- vs seqproc ann+edit (best config) ---")
        print(f"  seqproc ann+edit: {len(sp_ids):>10,} ({len(sp_ids)/1_000_000*100:.1f}%)")
        print(f"  splitcode dual:   {len(combined_ids):>10,} ({len(combined_ids)/1_000_000*100:.1f}%)")
        delta = len(sp_ids) - len(combined_ids)
        if delta > 0:
            print(f"  seqproc leads by: {delta:>10,} reads")
        else:
            print(f"  splitcode leads by: {-delta:>9,} reads")

    # Save results
    results = {
        "analysis": "splitcode_dual_orientation_supplementary",
        "dataset": "LR-SPLiT-seq (SRR13948564, 1M reads)",
        "note": "Simulated dual-orientation by running splitcode on forward + RC reads. "
                "Splitcode does not natively support orientation-aware matching. "
                "Combined runtime would be 2x a single run.",
        "forward_reads": len(fw_ids),
        "rc_reads": len(rc_ids),
        "both_orientations": len(both_ori),
        "forward_only": len(fw_only),
        "rc_only": len(rc_only),
        "combined_reads": len(combined_ids),
        "combined_recovery_pct": round(len(combined_ids) / 1_000_000 * 100, 2),
        "rc_runtime_s": round(rc_runtime, 2),
        "rc_memory_mb": round(rc_mem, 1),
    }

    if sp_ids:
        results["concordance_vs_seqproc"] = {
            "jaccard": round(jaccard(combined_ids, sp_ids), 4),
            "both": len(combined_ids & sp_ids),
            "splitcode_dual_only": len(combined_ids - sp_ids),
            "seqproc_only": len(sp_ids - combined_ids),
        }
    if mb_ids:
        results["concordance_vs_matchbox"] = {
            "jaccard": round(jaccard(combined_ids, mb_ids), 4),
            "both": len(combined_ids & mb_ids),
            "splitcode_dual_only": len(combined_ids - mb_ids),
            "matchbox_only": len(mb_ids - combined_ids),
        }

    out_json = RESULTS_DIR / "splitcode_dual_results.json"
    with open(out_json, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_json}")

    # Save combined IDs
    combined_ids_file = RESULTS_DIR / "splitcode_dual_ids.txt"
    with open(combined_ids_file, 'w') as f:
        for rid in sorted(combined_ids):
            f.write(rid + '\n')
    print(f"Saved: {combined_ids_file}")


if __name__ == "__main__":
    main()
