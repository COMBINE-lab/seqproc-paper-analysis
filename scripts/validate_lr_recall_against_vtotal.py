#!/usr/bin/env python3
"""
Compute |tool_emit_ids ∩ V_total_ids| / |V_total_ids| for each LR-SPLiT-seq tool.

This is the recall metric the seqproc paper headline reports: of the reads
that the independent structural-validity script accepts on the raw input
(V_total), what fraction did each tool recover?

  V_total      = reads in raw SRR13948564_full.fastq accepted by the strict
                 validity script (SplitSeqSingleEndValidityAnalyzer)
  numerator    = |tool's emitted read IDs  ∩  V_total IDs|
  denominator  = |V_total|
  Recall (%)   = 100 * numerator / denominator

For tools that are strictly more permissive than the validator
(edit-distance + orientation-aware), the intersection equals V_total and
recall is 100%. For forward-only tools, recall is roughly halved.

Output: JSON with exact counts and recall percentages for seqproc,
matchbox, splitcode (forward-only), and splitcode (dual-pass).

Inputs consumed (all on cluster, no re-benchmarking needed):
  - Raw FASTQ:           $SEQPROC_DATA_DIR/SRR13948564_full.fastq
  - seqproc emit IDs:    results/concordance/lr_splitseq/seqproc_edit_ids.txt
  - matchbox emit IDs:   results/concordance/lr_splitseq/matchbox_ids.txt
  - splitcode fw IDs:    results/concordance/lr_splitseq/splitcode_ids.txt
  - splitcode dual FQ:   $OUTDIR/splitcode_dual_combined_out.fq (from
                         splitcode_lr_dual_validate.py)

Usage:
    python3 scripts/validate_lr_recall_against_vtotal.py \\
        --dual-pass-fq results/splitcode_lr_dual/splitcode_dual_combined_out.fq \\
        --outdir results/lr_recall_vtotal

Expected runtime: ~3 minutes (V_total recompute ~150s + intersection ~30s).
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(
    os.environ.get("SEQPROC_PROJECT_ROOT", Path(__file__).resolve().parent.parent)
)
DATA_DIR = Path(os.environ.get("SEQPROC_DATA_DIR", PROJECT_ROOT / "data"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

BC1_WL = PROJECT_ROOT / "configs" / "seqproc" / "splitseq_bc1_seq2seq.tsv"
BC2_WL = PROJECT_ROOT / "configs" / "seqproc" / "splitseq_bc2_seq2seq.tsv"
BC3_WL = PROJECT_ROOT / "configs" / "seqproc" / "splitseq_bc3_seq2seq.tsv"

LR_INPUT = DATA_DIR / "SRR13948564_full.fastq"
CONCORDANCE_DIR = PROJECT_ROOT / "results" / "concordance" / "lr_splitseq"


def load_ids_txt(path: Path) -> set:
    """Load read IDs from a one-per-line text file."""
    if not path.exists():
        print(f"  [WARN] missing: {path}")
        return set()
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def fastq_ids(path: Path) -> set:
    """Extract read IDs from a FASTQ (first whitespace token, '@' stripped)."""
    ids = set()
    with open(path) as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                ids.add(line.strip().split()[0].lstrip("@"))
    return ids


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dual-pass-fq", required=True,
                        help="Path to splitcode_dual_combined_out.fq produced by "
                             "splitcode_lr_dual_validate.py")
    parser.add_argument("--outdir", required=True,
                        help="Output directory for V_total IDs and final JSON")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Disable validity caching so analyzer always recomputes (we need fresh IDs)
    import run_paper_benchmarks as rpb
    rpb._load_validity_cache = lambda *a, **kw: None
    rpb._save_validity_cache = lambda *a, **kw: None
    from run_paper_benchmarks import SplitSeqSingleEndValidityAnalyzer

    print("=" * 76)
    print("LR-SPLiT-seq RECALL vs V_total (independent strict validator)")
    print("=" * 76)
    print(f"  Raw FASTQ:        {LR_INPUT}")
    print(f"  Dual-pass FASTQ:  {args.dual_pass_fq}")
    print(f"  Output dir:       {outdir}")
    print()

    # ---- Step 1: V_total = validate raw input, save IDs ----
    print("[1/4] Running strict validator on raw input (computes V_total)...")
    analyzer = SplitSeqSingleEndValidityAnalyzer(
        str(BC1_WL), str(BC2_WL), str(BC3_WL)
    )
    t0 = time.time()
    v_total = analyzer.analyze_fastqs(str(LR_INPUT))
    print(f"  V_total = {len(v_total):,} reads  ({time.time()-t0:.1f}s)")

    # Save V_total IDs for future reuse
    v_total_path = outdir / "v_total_ids.txt"
    with open(v_total_path, "w") as f:
        for rid in sorted(v_total):
            f.write(rid + "\n")
    print(f"  Saved: {v_total_path}")

    # ---- Step 2: load per-tool emit ID sets ----
    print("\n[2/4] Loading per-tool emit ID sets from concordance cache...")
    seqproc_ids        = load_ids_txt(CONCORDANCE_DIR / "seqproc_edit_ids.txt")
    matchbox_ids       = load_ids_txt(CONCORDANCE_DIR / "matchbox_ids.txt")
    splitcode_fw_ids   = load_ids_txt(CONCORDANCE_DIR / "splitcode_ids.txt")
    print(f"  seqproc:        {len(seqproc_ids):>10,} emit IDs")
    print(f"  matchbox:       {len(matchbox_ids):>10,} emit IDs")
    print(f"  splitcode fw:   {len(splitcode_fw_ids):>10,} emit IDs")

    print("\n[3/4] Extracting splitcode dual-pass emit IDs from combined FASTQ...")
    t0 = time.time()
    splitcode_dual_ids = fastq_ids(Path(args.dual_pass_fq))
    print(f"  splitcode dual: {len(splitcode_dual_ids):>10,} emit IDs  "
          f"({time.time()-t0:.1f}s)")

    # ---- Step 3: compute intersections ----
    print("\n[4/4] Computing intersections with V_total...")
    tools = [
        ("seqproc",        seqproc_ids),
        ("matchbox",       matchbox_ids),
        ("splitcode_fwd",  splitcode_fw_ids),
        ("splitcode_dual", splitcode_dual_ids),
    ]
    out = {
        "v_total": len(v_total),
        "total_input_reads": 5_764_421,
        "validator": "SplitSeqSingleEndValidityAnalyzer "
                     "(both orientations, str.find linker, Hamming<=1 whitelist)",
        "tools": {},
    }
    print()
    print(f"{'Tool':<18} {'Emit':>12} {'∩ V_total':>12} {'Recall %':>10}  "
          f"{'V_total NOT in tool':>20}")
    print("-" * 76)
    for label, ids in tools:
        inter = v_total & ids
        missing = v_total - ids
        recall = 100 * len(inter) / len(v_total) if v_total else 0.0
        out["tools"][label] = {
            "emit": len(ids),
            "intersection_with_v_total": len(inter),
            "recall_pct": round(recall, 2),
            "v_total_reads_missed": len(missing),
        }
        print(f"{label:<18} {len(ids):>12,} {len(inter):>12,} {recall:>9.2f}%  "
              f"{len(missing):>20,}")

    out_json = outdir / "lr_recall_vtotal_results.json"
    with open(out_json, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {out_json}")
    print("\nThese are the EXACT numerator/denominator pairs for the LR-SPLiT-seq")
    print("row of Table 1 under the V_total-as-denominator framing:")
    print(f"  Denominator (V_total): {len(v_total):,}")
    print(f"  Numerator (tool ∩ V_total): see table above per tool.")


if __name__ == "__main__":
    main()
