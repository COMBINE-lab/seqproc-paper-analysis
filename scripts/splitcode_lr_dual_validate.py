#!/usr/bin/env python3
"""
Splitcode LR-SPLiT-seq dual-pass best-faith aggregate + validation.

This script answers two preprint-blocking questions for LR-SPLiT-seq splitcode:

  1. What is splitcode's recovery if we run it both forward and reverse-complement
     and take the union? (The existing table reports forward-only, which is
     unfair given the chemistry is randomly oriented.)
  2. What is splitcode's true per-output structural-validity precision on this
     chemistry? The published checkpoint reports 100% precision but a 1M-subset
     spot check showed ~17%, suggesting the validator may not have been run
     on splitcode's LR output in the original benchmark.

It does NOT re-run seqproc or matchbox. It does NOT re-run any other chemistry.
The only compute it consumes is two splitcode passes on LR-SPLiT-seq plus
~3 validation passes (output FASTQs + optional raw-input V_total).

Output: single JSON at $OUTDIR/splitcode_lr_dual_results.json with raw counts
for every cell needed to rebuild the LR-SPLiT-seq row of Table 1 honestly.

Usage on cluster (after sourcing setup_and_run.sh or equivalent):
    python3 scripts/splitcode_lr_dual_validate.py \\
        --dataset full \\
        --threads 8 \\
        --outdir $WORKDIR/results/splitcode_lr_dual

Env vars consumed (all set by setup_and_run.sh):
    SEQPROC_PROJECT_ROOT   path to this repo
    SEQPROC_DATA_DIR       directory holding SRR13948564_full.fastq
    SPLITCODE_BIN          absolute path to splitcode binary
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Set, Tuple

# ---------- path resolution (uses env vars set by setup_and_run.sh) ----------

PROJECT_ROOT = Path(
    os.environ.get("SEQPROC_PROJECT_ROOT", Path(__file__).resolve().parent.parent)
)
DATA_DIR = Path(
    os.environ.get("SEQPROC_DATA_DIR", PROJECT_ROOT / "data")
)
SPLITCODE_BIN = os.environ.get("SPLITCODE_BIN")
if not SPLITCODE_BIN:
    # Fall back to common locations
    for candidate in [
        PROJECT_ROOT.parent / "splitcode/build/src/splitcode",
        Path("/fs/nexus-projects/seqproc/bench/splitcode/build/src/splitcode"),
    ]:
        if candidate.exists():
            SPLITCODE_BIN = str(candidate)
            break
    if not SPLITCODE_BIN:
        sys.exit("ERROR: SPLITCODE_BIN env var not set and no fallback found")

SPLITCODE_CONFIG = PROJECT_ROOT / "configs" / "splitcode" / "splitseq_singleend.config"

BC1_WL = PROJECT_ROOT / "configs" / "seqproc" / "splitseq_bc1_seq2seq.tsv"
BC2_WL = PROJECT_ROOT / "configs" / "seqproc" / "splitseq_bc2_seq2seq.tsv"
BC3_WL = PROJECT_ROOT / "configs" / "seqproc" / "splitseq_bc3_seq2seq.tsv"

# Make the validator importable
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


# ---------- helpers ----------

def run_cmd(cmd: str) -> Tuple[float, float, int, str]:
    """Run shell cmd under /usr/bin/time -v. Returns (runtime_s, peak_mem_mb, exit, stderr)."""
    full_cmd = f"/usr/bin/time -v {cmd}"
    start = time.time()
    # Compat: capture_output/text were added in Python 3.7; spell them out for 3.6.
    result = subprocess.run(
        full_cmd, shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    runtime = time.time() - start
    peak_kb = 0
    exit_code = result.returncode
    for line in result.stderr.split("\n"):
        if "Maximum resident set size" in line:
            try:
                peak_kb = int(line.split(":")[1].strip())
            except (IndexError, ValueError):
                pass
        elif "Exit status" in line:
            try:
                exit_code = int(line.split(":")[1].strip())
            except (IndexError, ValueError):
                pass
    return runtime, peak_kb / 1024.0, exit_code, result.stderr


def count_fastq(path: Path) -> int:
    """Count FASTQ records (lines / 4)."""
    n = 0
    with open(path) as f:
        for i, _ in enumerate(f):
            if i % 4 == 0:
                n += 1
    return n


def read_ids(path: Path) -> Set[str]:
    """Extract read IDs (first whitespace-token, leading '@' stripped) from a FASTQ."""
    ids = set()
    with open(path) as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                ids.add(line.strip().split()[0].lstrip("@"))
    return ids


def rc_fastq_python(in_path: Path, out_path: Path) -> Tuple[int, float]:
    """Pure-Python reverse-complement (fallback). Returns (read_count, seconds)."""
    t0 = time.time()
    n = 0
    with open(in_path) as fin, open(out_path, "w") as fout:
        while True:
            h = fin.readline()
            if not h:
                break
            s = fin.readline().strip()
            p = fin.readline()
            q = fin.readline().strip()
            fout.write(h)
            fout.write(s.translate(COMPLEMENT)[::-1] + "\n")
            fout.write(p)
            fout.write(q[::-1] + "\n")
            n += 1
    return n, time.time() - t0


def rc_fastq(in_path: Path, out_path: Path, threads: int) -> Tuple[int, float, str]:
    """Try seqkit (fast, parallel) then fall back to Python. Returns (count, sec, method)."""
    if shutil.which("seqkit"):
        cmd = f"seqkit seq -r -p -t dna -j {threads} {in_path} -o {out_path}"
        runtime, _, exit_code, stderr = run_cmd(cmd)
        if exit_code == 0:
            return count_fastq(out_path), runtime, "seqkit"
        print(f"  seqkit failed (exit={exit_code}); falling back to python")
        print(stderr[-300:])
    n, t = rc_fastq_python(in_path, out_path)
    return n, t, "python"


def build_dedup_combined(fw_fq: Path, rc_fq: Path, out_fq: Path) -> Tuple[int, int]:
    """Write a deduplicated combined FASTQ (forward records, then RC-unique records).
    Returns (forward_count, total_combined_count)."""
    seen = set()
    fw_count = 0
    with open(out_fq, "w") as fout:
        # Forward pass — keep everything
        with open(fw_fq) as f:
            while True:
                h = f.readline()
                if not h:
                    break
                s = f.readline()
                p = f.readline()
                q = f.readline()
                rid = h.strip().split()[0].lstrip("@")
                seen.add(rid)
                fout.write(h + s + p + q)
                fw_count += 1
        # RC pass — only records not already in forward
        with open(rc_fq) as f:
            while True:
                h = f.readline()
                if not h:
                    break
                s = f.readline()
                p = f.readline()
                q = f.readline()
                rid = h.strip().split()[0].lstrip("@")
                if rid not in seen:
                    fout.write(h + s + p + q)
                    seen.add(rid)
    return fw_count, len(seen)


def run_splitcode(input_fq: Path, out_dir: Path, label: str, threads: int) -> dict:
    """Run splitcode and return {emitted, runtime_s, peak_mem_mb, out_fq}."""
    out_fq = out_dir / f"splitcode_{label}_out.fq"
    mapping = out_dir / f"splitcode_{label}_mapping.txt"
    cmd = (
        f"{SPLITCODE_BIN} -c {SPLITCODE_CONFIG} --assign -N 1 -t {threads} "
        f"-m {mapping} -o {out_fq} {input_fq}"
    )
    print(f"  Running splitcode ({label})... ", end="", flush=True)
    runtime, mem_mb, exit_code, stderr = run_cmd(cmd)
    if exit_code != 0:
        print(f"FAILED (exit={exit_code})")
        print(stderr[-500:])
        sys.exit(f"splitcode {label} pass failed")
    emitted = count_fastq(out_fq)
    print(f"emitted={emitted:,}  {runtime:.1f}s  {mem_mb:.0f}MB")
    return {
        "label": label,
        "emitted": emitted,
        "runtime_s": round(runtime, 2),
        "peak_mem_mb": round(mem_mb, 1),
        "out_fq": str(out_fq),
    }


def validate_fastq(path: Path, label: str) -> Tuple[int, float]:
    """Run SplitSeqSingleEndValidityAnalyzer; return (valid_count, seconds)."""
    # Import lazily so the validator class picks up monkeypatched cache dir
    from run_paper_benchmarks import SplitSeqSingleEndValidityAnalyzer
    import run_paper_benchmarks as rpb
    # Disable on-disk caching so each call recomputes against current FASTQ
    rpb._load_validity_cache = lambda *a, **kw: None
    rpb._save_validity_cache = lambda *a, **kw: None

    analyzer = SplitSeqSingleEndValidityAnalyzer(str(BC1_WL), str(BC2_WL), str(BC3_WL))
    print(f"  Validating {label}... ", end="", flush=True)
    t0 = time.time()
    valid = analyzer.analyze_fastqs(str(path))
    elapsed = time.time() - t0
    print(f"valid={len(valid):,}  ({elapsed:.1f}s)")
    return len(valid), elapsed


# ---------- main ----------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dataset", choices=["1m", "full"], default="full",
                        help="LR-SPLiT-seq subset to run on (default: full)")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--outdir", required=True,
                        help="Output directory for intermediate FASTQs and JSON")
    parser.add_argument("--skip-vtotal", action="store_true",
                        help="Skip raw-input validation (V_total). Saves ~3-5 min.")
    parser.add_argument("--reuse-rc", action="store_true",
                        help="Reuse existing RC FASTQ if it's already in outdir")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fastq_name = (
        "SRR13948564_1M.fastq" if args.dataset == "1m" else "SRR13948564_full.fastq"
    )
    fw_input = DATA_DIR / fastq_name
    if not fw_input.exists():
        sys.exit(f"ERROR: input FASTQ not found: {fw_input}")

    print("=" * 78)
    print("SPLITCODE LR-SPLiT-seq DUAL-PASS BEST-FAITH AGGREGATE + VALIDATION")
    print("=" * 78)
    print(f"  PROJECT_ROOT:    {PROJECT_ROOT}")
    print(f"  DATA_DIR:        {DATA_DIR}")
    print(f"  SPLITCODE_BIN:   {SPLITCODE_BIN}")
    print(f"  Input FASTQ:     {fw_input}")
    print(f"  Output dir:      {outdir}")
    print(f"  Threads:         {args.threads}")
    print(f"  Skip V_total:    {args.skip_vtotal}")
    print()

    # ---- step 1: count input reads ----
    print("[1/6] Counting input reads...")
    t0 = time.time()
    total_input = count_fastq(fw_input)
    print(f"  Total input reads: {total_input:,}  ({time.time()-t0:.1f}s)")

    # ---- step 2: forward splitcode ----
    print("\n[2/6] Splitcode FORWARD pass...")
    fw_result = run_splitcode(fw_input, outdir, "fw", args.threads)

    # ---- step 3: reverse-complement input ----
    print("\n[3/6] Reverse-complementing input FASTQ...")
    rc_input = outdir / f"{fastq_name.replace('.fastq', '_RC.fastq')}"
    if args.reuse_rc and rc_input.exists():
        print(f"  Reusing existing RC FASTQ: {rc_input}")
        rc_build_time = 0.0
        rc_method = "cached"
    else:
        n_rc, rc_build_time, rc_method = rc_fastq(fw_input, rc_input, args.threads)
        print(f"  Wrote {n_rc:,} records via {rc_method} in {rc_build_time:.1f}s")

    # ---- step 4: RC splitcode ----
    print("\n[4/6] Splitcode RC pass...")
    rc_result = run_splitcode(rc_input, outdir, "rc", args.threads)

    # ---- step 5: build dedup combined FASTQ + ID-set stats ----
    print("\n[5/6] Building deduplicated combined FASTQ...")
    fw_ids = read_ids(Path(fw_result["out_fq"]))
    rc_ids = read_ids(Path(rc_result["out_fq"]))
    both = fw_ids & rc_ids
    fw_only = fw_ids - rc_ids
    rc_only = rc_ids - fw_ids
    union = fw_ids | rc_ids
    combined_fq = outdir / "splitcode_dual_combined_out.fq"
    _, combined_count = build_dedup_combined(
        Path(fw_result["out_fq"]), Path(rc_result["out_fq"]), combined_fq
    )
    print(f"  Forward emit:   {len(fw_ids):>12,}")
    print(f"  RC emit:        {len(rc_ids):>12,}")
    print(f"  Both:           {len(both):>12,}")
    print(f"  Forward-unique: {len(fw_only):>12,}")
    print(f"  RC-unique:      {len(rc_only):>12,}")
    print(f"  Union (dual):   {len(union):>12,}  ({100*len(union)/total_input:.2f}% of input)")

    # ---- step 6: validation ----
    print("\n[6/6] Validation (independent structural-validity script)...")
    fw_valid, fw_val_t = validate_fastq(Path(fw_result["out_fq"]), "forward output")
    rc_valid, rc_val_t = validate_fastq(Path(rc_result["out_fq"]), "RC output")
    combined_valid, combined_val_t = validate_fastq(combined_fq, "dual-pass combined")
    if args.skip_vtotal:
        vtotal_valid, vtotal_val_t = None, 0.0
        print("  Skipping V_total (raw-input validation)")
    else:
        vtotal_valid, vtotal_val_t = validate_fastq(fw_input, "raw input (V_total)")

    summed_runtime = fw_result["runtime_s"] + rc_result["runtime_s"]

    print()
    print("=" * 78)
    print("RESULTS SUMMARY")
    print("=" * 78)
    print(f"Total input reads:                 {total_input:>12,}")
    print(f"Forward emit / valid / precision:  {len(fw_ids):>12,} / {fw_valid:>10,} / "
          f"{100*fw_valid/max(len(fw_ids),1):.1f}%")
    print(f"RC emit / valid / precision:       {len(rc_ids):>12,} / {rc_valid:>10,} / "
          f"{100*rc_valid/max(len(rc_ids),1):.1f}%")
    print(f"Dual emit / valid / precision:     {combined_count:>12,} / {combined_valid:>10,} / "
          f"{100*combined_valid/max(combined_count,1):.1f}%")
    if vtotal_valid is not None:
        print(f"V_total (raw-input valid):         {vtotal_valid:>12,}  "
              f"({100*vtotal_valid/total_input:.2f}% of input)")
    print(f"Summed splitcode runtime:          {summed_runtime:.1f}s "
          f"(fw {fw_result['runtime_s']:.1f}s + rc {rc_result['runtime_s']:.1f}s)")
    print(f"Peak memory (single pass):         {max(fw_result['peak_mem_mb'], rc_result['peak_mem_mb']):.0f} MB")
    print()
    print(f"OLD METRIC (valid / total_input):")
    print(f"  Forward-only:  {100*fw_valid/total_input:.2f}%")
    print(f"  Dual-pass:     {100*combined_valid/total_input:.2f}%")
    if vtotal_valid:
        print(f"NEW METRIC (valid / V_total):")
        print(f"  Forward-only:  {100*fw_valid/vtotal_valid:.2f}%")
        print(f"  Dual-pass:     {100*combined_valid/vtotal_valid:.2f}%")

    # ---- save JSON ----
    results = {
        "analysis": "splitcode_lr_dual_validate",
        "dataset": fastq_name,
        "input_path": str(fw_input),
        "total_input_reads": total_input,
        "splitcode_binary": SPLITCODE_BIN,
        "splitcode_config": str(SPLITCODE_CONFIG),
        "validator": "SplitSeqSingleEndValidityAnalyzer (Hamming <=1 on barcodes, exact str.find for linker)",
        "threads": args.threads,
        "forward": {
            "emitted": len(fw_ids),
            "valid": fw_valid,
            "precision_pct": round(100 * fw_valid / max(len(fw_ids), 1), 2),
            "runtime_s": fw_result["runtime_s"],
            "peak_mem_mb": fw_result["peak_mem_mb"],
        },
        "rc": {
            "emitted": len(rc_ids),
            "valid": rc_valid,
            "precision_pct": round(100 * rc_valid / max(len(rc_ids), 1), 2),
            "runtime_s": rc_result["runtime_s"],
            "peak_mem_mb": rc_result["peak_mem_mb"],
            "rc_build_seconds": round(rc_build_time, 2),
            "rc_build_method": rc_method,
        },
        "dual_pass": {
            "emitted_union": combined_count,
            "valid": combined_valid,
            "precision_pct": round(100 * combined_valid / max(combined_count, 1), 2),
            "both_orientations": len(both),
            "forward_unique": len(fw_only),
            "rc_unique": len(rc_only),
            "summed_runtime_s": round(summed_runtime, 2),
            "peak_mem_mb_single_pass": max(fw_result["peak_mem_mb"], rc_result["peak_mem_mb"]),
        },
        "v_total": {
            "valid": vtotal_valid,
            "pct_of_input": round(100 * vtotal_valid / total_input, 2) if vtotal_valid else None,
        } if vtotal_valid is not None else None,
        "old_metric_valid_over_input": {
            "forward_only": round(100 * fw_valid / total_input, 2),
            "dual_pass": round(100 * combined_valid / total_input, 2),
        },
        "new_metric_valid_over_vtotal": (
            {
                "forward_only": round(100 * fw_valid / vtotal_valid, 2),
                "dual_pass": round(100 * combined_valid / vtotal_valid, 2),
            }
            if vtotal_valid else None
        ),
    }

    out_json = outdir / "splitcode_lr_dual_results.json"
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_json}")


if __name__ == "__main__":
    main()
