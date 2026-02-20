#!/usr/bin/env python3
"""
Phase 5 Pre-requisite: Re-run LR-SPLiT-seq performance benchmarks.

Sprint 4 used a forward-only seqproc config for LR-SPLiT-seq (23.8% recovery,
~2.1s). The annotation+edit config (49.9% recovery, ~5.1s single run) is the
correct config for the paper. This script runs 3 replicates for all tools
to get proper mean/std for Table 2.

Usage:
    python3 scripts/phase5_lr_perf_rerun.py --threads 4 --reps 3
"""

import json
import os
import subprocess
import time
import argparse
import statistics
from pathlib import Path
from typing import Tuple

PROJECT_ROOT = Path(__file__).parent.parent
CONFIGS = PROJECT_ROOT / "configs"
RESULTS_DIR = PROJECT_ROOT / "results" / "phase5_lr_perf"

SEQPROC_BIN = os.environ.get(
    "SEQPROC_BIN",
    str(PROJECT_ROOT.parent / "combine-lab/seqproc/target/release/seqproc")
)
MATCHBOX_BIN = os.environ.get(
    "MATCHBOX_BIN",
    str(PROJECT_ROOT.parent / "matchbox/target/release/matchbox")
)
SPLITCODE_BIN = os.environ.get(
    "SPLITCODE_BIN",
    str(PROJECT_ROOT.parent / "splitcode/build/src/splitcode")
)

LR_FASTQ = PROJECT_ROOT / "data" / "SRR13948564_1M.fastq"


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


def count_fastq_reads(filepath: str) -> int:
    """Count reads in a FASTQ file."""
    count = 0
    if not os.path.exists(filepath):
        return 0
    with open(filepath) as f:
        for _ in f:
            count += 1
    return count // 4


def main():
    parser = argparse.ArgumentParser(description='LR-SPLiT-seq 3-rep performance benchmark')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--reps', type=int, default=3)
    args = parser.parse_args()

    print("=" * 70)
    print("PHASE 5 PRE-REQ: LR-SPLiT-seq Performance Re-run")
    print(f"  Replicates: {args.reps}, Threads: {args.threads}")
    print("=" * 70)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    if not LR_FASTQ.exists():
        print(f"[ERROR] Input FASTQ not found: {LR_FASTQ}")
        return

    results = {"dataset": "LR-SPLiT-seq (SRR13948564, 1M reads)",
               "threads": args.threads, "replicates": args.reps, "tools": {}}

    # ---- seqproc (annotation + edit distance) ----
    print("\n--- seqproc (ann+edit) ---")
    geom = CONFIGS / "seqproc/splitseq_singleend_edit_ann.geom"

    sp_runs = []
    for rep in range(1, args.reps + 1):
        out_fq = RESULTS_DIR / f"seqproc_rep{rep}_R1.fq"
        cmd = (f"{SEQPROC_BIN} --geom {geom} --file1 {LR_FASTQ} "
               f"--threads {args.threads} --out1 {out_fq}")
        runtime, mem, rc, stderr = run_cmd(cmd, PROJECT_ROOT)
        reads_out = count_fastq_reads(str(out_fq))
        print(f"  Rep {rep}: {runtime:.2f}s, {mem:.1f}MB, {reads_out:,} reads, rc={rc}")
        if rc != 0:
            print(f"    [ERROR] {stderr[-300:]}")
        sp_runs.append({"runtime": runtime, "memory_mb": mem, "reads_out": reads_out, "rc": rc})

    runtimes = [r["runtime"] for r in sp_runs if r["rc"] == 0]
    mems = [r["memory_mb"] for r in sp_runs if r["rc"] == 0]
    reads = [r["reads_out"] for r in sp_runs if r["rc"] == 0]
    results["tools"]["seqproc"] = {
        "config": "splitseq_singleend_edit_ann.geom (ann+edit)",
        "runs": sp_runs,
        "mean_runtime": statistics.mean(runtimes) if runtimes else 0,
        "std_runtime": statistics.stdev(runtimes) if len(runtimes) > 1 else 0,
        "mean_memory_mb": statistics.mean(mems) if mems else 0,
        "mean_reads_out": statistics.mean(reads) if reads else 0,
        "recovery_pct": round(statistics.mean(reads) / 1_000_000 * 100, 2) if reads else 0,
    }
    print(f"  Mean: {results['tools']['seqproc']['mean_runtime']:.2f}s "
          f"+/- {results['tools']['seqproc']['std_runtime']:.2f}s, "
          f"{results['tools']['seqproc']['mean_memory_mb']:.1f}MB, "
          f"{results['tools']['seqproc']['recovery_pct']}% recovery")

    # ---- matchbox (dual-orientation) ----
    print("\n--- matchbox (dual) ---")
    mb_config = CONFIGS / "matchbox/splitseq_singleend_dual.mb"
    mb_runs = []
    for rep in range(1, args.reps + 1):
        out_tsv = RESULTS_DIR / f"matchbox_rep{rep}.tsv"
        cmd = (f"{MATCHBOX_BIN} -e 0.2 -t {args.threads} "
               f"-s {mb_config} {LR_FASTQ} > {out_tsv}")
        runtime, mem, rc, stderr = run_cmd(cmd, PROJECT_ROOT)
        # Count matchbox output from TSV (no header, 1 line = 1 read)
        reads_out = 0
        if out_tsv.exists():
            with open(out_tsv) as f:
                reads_out = sum(1 for _ in f)
        # Clean any FASTQ files matchbox may have written
        for fq in ['mb_r1.fq', 'mb_r2.fq']:
            src = PROJECT_ROOT / fq
            if src.exists():
                src.unlink()
        print(f"  Rep {rep}: {runtime:.2f}s, {mem:.1f}MB, {reads_out:,} reads, rc={rc}")
        if rc != 0:
            print(f"    [ERROR] {stderr[-300:]}")
        mb_runs.append({"runtime": runtime, "memory_mb": mem, "reads_out": reads_out, "rc": rc})

    runtimes = [r["runtime"] for r in mb_runs if r["rc"] == 0]
    mems = [r["memory_mb"] for r in mb_runs if r["rc"] == 0]
    reads = [r["reads_out"] for r in mb_runs if r["rc"] == 0]
    results["tools"]["matchbox"] = {
        "config": "splitseq_singleend_dual.mb",
        "runs": mb_runs,
        "mean_runtime": statistics.mean(runtimes) if runtimes else 0,
        "std_runtime": statistics.stdev(runtimes) if len(runtimes) > 1 else 0,
        "mean_memory_mb": statistics.mean(mems) if mems else 0,
        "mean_reads_out": statistics.mean(reads) if reads else 0,
        "recovery_pct": round(statistics.mean(reads) / 1_000_000 * 100, 2) if reads else 0,
    }
    print(f"  Mean: {results['tools']['matchbox']['mean_runtime']:.2f}s "
          f"+/- {results['tools']['matchbox']['std_runtime']:.2f}s, "
          f"{results['tools']['matchbox']['mean_memory_mb']:.1f}MB, "
          f"{results['tools']['matchbox']['recovery_pct']}% recovery")

    # ---- splitcode (forward-only, no native orientation support) ----
    print("\n--- splitcode (forward-only) ---")
    sc_config = CONFIGS / "splitcode/splitseq_singleend.config"
    sc_runs = []
    for rep in range(1, args.reps + 1):
        sc_out = RESULTS_DIR / f"splitcode_rep{rep}.fq"
        sc_map = RESULTS_DIR / f"splitcode_rep{rep}_mapping.txt"
        cmd = (f"{SPLITCODE_BIN} -c {sc_config} "
               f"--assign -N 1 -t {args.threads} -m {sc_map} "
               f"-o {sc_out} {LR_FASTQ}")
        runtime, mem, rc, stderr = run_cmd(cmd, PROJECT_ROOT)
        reads_out = count_fastq_reads(str(sc_out))
        print(f"  Rep {rep}: {runtime:.2f}s, {mem:.1f}MB, {reads_out:,} reads, rc={rc}")
        if rc != 0:
            print(f"    [ERROR] {stderr[-300:]}")
        sc_runs.append({"runtime": runtime, "memory_mb": mem, "reads_out": reads_out, "rc": rc})

    runtimes = [r["runtime"] for r in sc_runs if r["rc"] == 0]
    mems = [r["memory_mb"] for r in sc_runs if r["rc"] == 0]
    reads = [r["reads_out"] for r in sc_runs if r["rc"] == 0]
    results["tools"]["splitcode"] = {
        "config": "splitseq_singleend.config (forward-only)",
        "runs": sc_runs,
        "mean_runtime": statistics.mean(runtimes) if runtimes else 0,
        "std_runtime": statistics.stdev(runtimes) if len(runtimes) > 1 else 0,
        "mean_memory_mb": statistics.mean(mems) if mems else 0,
        "mean_reads_out": statistics.mean(reads) if reads else 0,
        "recovery_pct": round(statistics.mean(reads) / 1_000_000 * 100, 2) if reads else 0,
    }
    print(f"  Mean: {results['tools']['splitcode']['mean_runtime']:.2f}s "
          f"+/- {results['tools']['splitcode']['std_runtime']:.2f}s, "
          f"{results['tools']['splitcode']['mean_memory_mb']:.1f}MB, "
          f"{results['tools']['splitcode']['recovery_pct']}% recovery")

    # ---- Summary ----
    print(f"\n{'='*70}")
    print("SUMMARY: LR-SPLiT-seq Performance (1M reads, annotation configs)")
    print(f"{'='*70}")
    print(f"{'Tool':<15} {'Runtime':>15} {'Memory':>12} {'Recovery':>10}")
    print("-" * 55)
    for tool in ['seqproc', 'matchbox', 'splitcode']:
        t = results['tools'][tool]
        rt = f"{t['mean_runtime']:.1f} +/- {t['std_runtime']:.1f}s"
        mem = f"{t['mean_memory_mb']:.0f}MB"
        rec = f"{t['recovery_pct']}%"
        print(f"{tool:<15} {rt:>15} {mem:>12} {rec:>10}")

    # Save results
    out_json = RESULTS_DIR / "lr_splitseq_perf_results.json"
    with open(out_json, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_json}")


if __name__ == "__main__":
    main()
