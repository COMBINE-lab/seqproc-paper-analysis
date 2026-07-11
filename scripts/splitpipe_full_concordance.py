#!/usr/bin/env python3
"""
Full-dataset concordance of seqproc / matchbox / splitcode against the vendor
tool split-pipe, on SPLiT-seq paired-end data (SRR6750041, full 86.8M reads).

This computes the same precision / recall / F1 / Jaccard that the paper already
reports for seqproc-vs-split-pipe on the 10M subset, but on the FULL dataset and
for ALL THREE tools -- in particular giving splitcode's precision against the
vendor reference (the number the paper currently only infers).

Ground truth = the set of reads split-pipe accepts as having a valid barcode.
split-pipe writes these to <run>/process/barcode_head.fastq, where each header
embeds the original read ID after the final "@", e.g.
    @04_33_..._OH_@SRR6750041.71 1   ->   SRR6750041.71
That ID space matches the tool emit-ID files in
results/concordance/splitseq_pe/{seqproc,matchbox,splitcode}_ids.txt
(extracted as header.split()[0].lstrip('@')).

For each tool T:
    inter     = |T_emit_ids  ∩  splitpipe_valid_ids|
    precision = inter / |T_emit_ids|       (fraction of T's output that is valid)
    recall    = inter / |splitpipe_valid|  (fraction of vendor-valid reads T found)
    F1        = 2PR / (P + R)
    Jaccard   = inter / |T_emit_ids ∪ splitpipe_valid_ids|

PREREQUISITE (the gating step -- run on the cluster first):
    split-pipe must be run on the FULL SRR6750041 R1/R2 (header-fixed), with the
    same chemistry/kit as the existing 10M run. Only the barcode-detection output
    (barcode_head.fastq) is needed; the genome/alignment is irrelevant here.

Usage:
    python3 scripts/splitpipe_full_concordance.py \\
        --splitpipe-fastq /path/to/splitpipe_results_FULL/process/barcode_head.fastq \\
        --concordance-dir results/concordance/splitseq_pe \\
        --outdir results/splitpipe_full_concordance
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path


def extract_splitpipe_ids(fastq_path: Path) -> set:
    """Read split-pipe barcode_head.fastq; return set of original read IDs.

    Each header is '@<barcode-prefix>_@<ORIG_ID> <rest>'. The original read ID is
    the token after the final '@', up to the first whitespace.
    """
    ids = set()
    n = 0
    with open(fastq_path) as f:
        for i, line in enumerate(f):
            if i % 4 != 0:
                continue
            n += 1
            if n % 5_000_000 == 0:
                print(f"    ...{n:,} split-pipe records parsed, {len(ids):,} ids")
            at = line.rfind("@")
            if at < 0:
                continue
            ids.add(line[at + 1:].split()[0])
    return ids


def load_tool_ids(path: Path) -> set:
    if not path.exists():
        return set()
    with open(path) as f:
        return {ln.strip() for ln in f if ln.strip()}


def find_tool_file(concordance_dir: Path, names) -> Path:
    for nm in names:
        p = concordance_dir / nm
        if p.exists():
            return p
    return concordance_dir / names[0]  # report the expected name even if missing


def metrics(tool_ids: set, gt_ids: set) -> dict:
    inter = len(tool_ids & gt_ids)
    union = len(tool_ids | gt_ids)
    emit = len(tool_ids)
    gt = len(gt_ids)
    precision = inter / emit if emit else 0.0
    recall = inter / gt if gt else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0
    jaccard = inter / union if union else 0.0
    return {
        "emitted": emit,
        "intersection_with_splitpipe": inter,
        "precision": round(100 * precision, 2),
        "recall": round(100 * recall, 2),
        "f1": round(f1, 4),
        "jaccard": round(jaccard, 4),
    }


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--splitpipe-fastq", required=True,
                    help="barcode_head.fastq from the FULL-dataset split-pipe run")
    ap.add_argument("--concordance-dir",
                    default="results/concordance/splitseq_pe",
                    help="dir holding {seqproc,matchbox,splitcode}_ids.txt")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    sp_fastq = Path(args.splitpipe_fastq)
    cdir = Path(args.concordance_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not sp_fastq.exists():
        sys.exit(f"ERROR: split-pipe fastq not found: {sp_fastq}\n"
                 "Run split-pipe on the full dataset first (see module docstring).")

    print("=" * 76)
    print("FULL-DATASET split-pipe concordance (SPLiT-seq PE)")
    print("=" * 76)

    print(f"[1/3] Parsing split-pipe valid reads from {sp_fastq} ...")
    t0 = time.time()
    gt = extract_splitpipe_ids(sp_fastq)
    print(f"  split-pipe valid (ground truth): {len(gt):,}  ({time.time()-t0:.1f}s)")
    # cache the GT ids
    with open(outdir / "splitpipe_valid_ids.txt", "w") as fh:
        for rid in sorted(gt):
            fh.write(rid + "\n")

    print("\n[2/3] Loading tool emit-ID sets ...")
    tool_files = {
        "seqproc":   find_tool_file(cdir, ["seqproc_edit_ids.txt", "seqproc_ids.txt"]),
        "matchbox":  find_tool_file(cdir, ["matchbox_ids.txt"]),
        "splitcode": find_tool_file(cdir, ["splitcode_ids.txt"]),
    }
    tool_ids = {}
    for tool, p in tool_files.items():
        ids = load_tool_ids(p)
        tool_ids[tool] = ids
        flag = "" if ids else "  [MISSING -- run concordance_analysis.py first]"
        print(f"  {tool:<10} {len(ids):>12,} emit ids  ({p}){flag}")

    print("\n[3/3] Computing metrics vs split-pipe ...\n")
    out = {
        "dataset": "SRR6750041 (SPLiT-seq PE, full)",
        "splitpipe_valid": len(gt),
        "splitpipe_fastq": str(sp_fastq),
        "tools": {},
    }
    print(f"{'Tool':<10} {'Emit':>13} {'∩ split-pipe':>14} {'Prec%':>7} "
          f"{'Recall%':>8} {'F1':>7} {'Jaccard':>8}")
    print("-" * 76)
    for tool, ids in tool_ids.items():
        if not ids:
            continue
        m = metrics(ids, gt)
        out["tools"][tool] = m
        print(f"{tool:<10} {m['emitted']:>13,} {m['intersection_with_splitpipe']:>14,} "
              f"{m['precision']:>7.2f} {m['recall']:>8.2f} {m['f1']:>7.4f} "
              f"{m['jaccard']:>8.4f}")

    out_json = outdir / "splitpipe_full_concordance.json"
    with open(out_json, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nSaved: {out_json}")
    print("\nPrecision = fraction of each tool's output that split-pipe also calls valid.")
    print("Recall    = fraction of split-pipe's valid reads each tool recovered.")


if __name__ == "__main__":
    main()
