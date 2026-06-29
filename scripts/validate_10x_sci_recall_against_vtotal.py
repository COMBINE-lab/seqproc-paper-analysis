#!/usr/bin/env python3
"""
Compute recall against V_total for 10x Chromium v2 and sci-RNA-seq3, so the
Recovery column can use one denominator (recall against V_total) across the
three well-behaved chemistries (PE, 10x, sci). LR-SPLiT-seq is handled
separately because its strict exact-match validator is not comparable to the
tools' edit-distance matching.

For the chosen chemistry it:
  1. runs the appropriate structural-validity analyzer on the RAW input to get
     V_total (TenXValidityAnalyzer for 10x, SciSeqValidityAnalyzer for sci),
  2. loads each tool's emit-ID set from results/concordance/<chem>/,
  3. reports |tool emit ∩ V_total| / |V_total| (recall) per tool.

Expectation: 10x lands at ~100% for all tools (trivial fixed geometry);
sci lands near the existing validated-emitted figure (~99%+).

Usage (cluster, full data, after concordance_analysis.py has populated IDs):
    python3 scripts/validate_10x_sci_recall_against_vtotal.py --chemistry 10x \
        --outdir results/vtotal_10x
    python3 scripts/validate_10x_sci_recall_against_vtotal.py --chemistry sci \
        --outdir results/vtotal_sci

Env vars (set by setup_and_run.sh): SEQPROC_PROJECT_ROOT, SEQPROC_DATA_DIR.
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(os.environ.get("SEQPROC_PROJECT_ROOT", Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

import data_config  # noqa: E402
import run_paper_benchmarks as rpb  # noqa: E402

# always recompute, do not read a stale cache
rpb._load_validity_cache = lambda *a, **kw: None
rpb._save_validity_cache = lambda *a, **kw: None


def load_ids(path: Path) -> set:
    if not path.exists():
        return set()
    with open(path) as f:
        return {ln.strip() for ln in f if ln.strip()}


def find_tool_file(cdir: Path, names) -> Path:
    for nm in names:
        if (cdir / nm).exists():
            return cdir / nm
    return cdir / names[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chemistry", required=True, choices=["10x", "sci"])
    ap.add_argument("--reads", default="full", choices=["1m", "full"])
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    datasets = data_config.resolve_datasets(args.reads)
    if args.chemistry == "10x":
        ds_key, cdir_name = "10x_short", "10x_short"
        r1, r2, total = datasets["10x_short"][args.reads] if args.reads in datasets["10x_short"] else (None, None, None)
        analyzer = rpb.TenXValidityAnalyzer(is_short_read=True)
        gt = lambda: analyzer.analyze_fastqs(str(r1), str(r2))
    else:
        ds_key, cdir_name = "sciseq", "sciseq"
        r1, r2, total = datasets["sciseq"][args.reads] if args.reads in datasets["sciseq"] else (None, None, None)
        analyzer = rpb.SciSeqValidityAnalyzer()
        gt = lambda: analyzer.analyze_fastqs(str(r1), str(r2))

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cdir = PROJECT_ROOT / "results" / "concordance" / cdir_name

    print("=" * 74)
    print(f"V_total recall  --  {args.chemistry}  ({args.reads})")
    print("=" * 74)
    print(f"  R1: {r1}\n  R2: {r2}\n  concordance dir: {cdir}\n")

    print("[1/3] Running validator on raw input to build V_total ...")
    t0 = time.time()
    v_total = set(gt())
    print(f"  V_total = {len(v_total):,}  ({time.time()-t0:.1f}s)")
    with open(outdir / "v_total_ids.txt", "w") as fh:
        for rid in sorted(v_total):
            fh.write(rid + "\n")

    print("\n[2/3] Loading tool emit-ID sets ...")
    tool_files = {
        "seqproc":  find_tool_file(cdir, ["seqproc_edit_ids.txt", "seqproc_ids.txt"]),
        "matchbox": find_tool_file(cdir, ["matchbox_ids.txt"]),
        "splitcode": find_tool_file(cdir, ["splitcode_ids.txt"]),
    }

    print("\n[3/3] Recall against V_total ...\n")
    out = {"chemistry": args.chemistry, "v_total": len(v_total), "tools": {}}
    print(f"{'Tool':<11}{'Emit':>12}{'∩ V_total':>12}{'Recall%':>9}")
    print("-" * 44)
    for tool, p in tool_files.items():
        ids = load_ids(p)
        if not ids:
            print(f"{tool:<11}{'MISSING':>12}   run concordance_analysis.py first")
            continue
        inter = len(ids & v_total)
        recall = 100 * inter / len(v_total) if v_total else 0
        out["tools"][tool] = {"emit": len(ids), "intersection": inter, "recall_pct": round(recall, 2)}
        print(f"{tool:<11}{len(ids):>12,}{inter:>12,}{recall:>8.2f}%")

    with open(outdir / f"vtotal_recall_{args.chemistry}.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nSaved: {outdir / f'vtotal_recall_{args.chemistry}.json'}")


if __name__ == "__main__":
    main()
