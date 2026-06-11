#!/usr/bin/env python3
"""
Compute |tool_emit_ids ∩ V_total_ids| / |V_total_ids| for each SPLiT-seq PE tool.

This is the same metric the paper's LR-SPLiT-seq row of Table 3 already uses,
extended to SPLiT-seq PE for cross-row consistency.

  V_total      = reads in raw SRR6750041_R2.fastq (where the linker lives)
                 accepted by the strict validity script
                 (SplitSeqValidityAnalyzer)
  numerator    = |tool's emitted read IDs ∩ V_total IDs|
  denominator  = |V_total|
  Recall (%)   = 100 * numerator / denominator

Expected outcome: seqproc and splitcode (which both preserve linkers and
use comparable matching criteria) should sit near 100%. matchbox may sit
slightly lower if its edit-distance grammar diverges from the validator
on a meaningful slice (same pattern as LR-SPLiT-seq, where matchbox
came in at 87.6%).

Inputs (all on cluster, no benchmark re-run):
  - Raw R1 FASTQ: $SEQPROC_DATA_DIR/SRR6750041_R1.fastq
  - Raw R2 FASTQ: $SEQPROC_DATA_DIR/SRR6750041_R2.fastq
  - seqproc emit IDs:  results/concordance/splitseq_pe/seqproc_ids.txt
                       (or seqproc_edit_ids.txt -- the script tries both)
  - matchbox emit IDs: results/concordance/splitseq_pe/matchbox_ids.txt
  - splitcode emit IDs: results/concordance/splitseq_pe/splitcode_ids.txt

Output: JSON at $OUTDIR/pe_recall_vtotal_results.json with exact counts
and recall percentages.

Usage:
    python3 scripts/validate_pe_recall_against_vtotal.py \\
        --outdir /fs/nexus-projects/seqproc/bench/results/pe_recall_vtotal

Expected runtime: ~25-40 min (V_total recompute on 86.8M PE reads is the
bottleneck; intersections are ~30 seconds).
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

PE_R1 = DATA_DIR / "SRR6750041_R1.fastq"
PE_R2 = DATA_DIR / "SRR6750041_R2.fastq"
CONCORDANCE_DIR = PROJECT_ROOT / "results" / "concordance" / "splitseq_pe"


def load_ids_txt(path: Path) -> set:
    """Load read IDs from a one-per-line text file."""
    if not path.exists():
        return set()
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def find_seqproc_ids() -> tuple:
    """seqproc PE may use 'seqproc_ids.txt' or 'seqproc_edit_ids.txt'. Try both."""
    for name in ("seqproc_edit_ids.txt", "seqproc_ids.txt"):
        p = CONCORDANCE_DIR / name
        if p.exists():
            return name, load_ids_txt(p)
    return None, set()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", required=True,
                        help="Output directory for V_total IDs and final JSON")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Disable validity caching to force a fresh recompute
    import run_paper_benchmarks as rpb
    rpb._load_validity_cache = lambda *a, **kw: None
    rpb._save_validity_cache = lambda *a, **kw: None
    from run_paper_benchmarks import SplitSeqValidityAnalyzer

    print("=" * 76)
    print("SPLiT-seq PE RECALL vs V_total (independent strict validator)")
    print("=" * 76)
    print(f"  Raw R1:   {PE_R1}")
    print(f"  Raw R2:   {PE_R2}")
    print(f"  Concordance ID cache: {CONCORDANCE_DIR}")
    print(f"  Output dir: {outdir}")
    print()

    for p in (PE_R1, PE_R2):
        if not p.exists():
            sys.exit(f"ERROR: missing input FASTQ: {p}")

    # ---- Step 1: compute V_total on raw input ----
    print("[1/4] Running strict validator on raw input (computes V_total)...")
    analyzer = SplitSeqValidityAnalyzer(
        str(BC1_WL), str(BC2_WL), str(BC3_WL)
    )
    t0 = time.time()
    v_total = analyzer.analyze_fastqs(str(PE_R1), str(PE_R2))
    print(f"  V_total = {len(v_total):,} reads  ({time.time()-t0:.1f}s)")

    v_total_path = outdir / "v_total_pe_ids.txt"
    with open(v_total_path, "w") as f:
        for rid in sorted(v_total):
            f.write(rid + "\n")
    print(f"  Saved: {v_total_path}")

    # ---- Step 2: load per-tool emit ID sets ----
    print("\n[2/4] Loading per-tool emit ID sets...")
    seqproc_name, seqproc_ids = find_seqproc_ids()
    matchbox_ids  = load_ids_txt(CONCORDANCE_DIR / "matchbox_ids.txt")
    splitcode_ids = load_ids_txt(CONCORDANCE_DIR / "splitcode_ids.txt")
    print(f"  seqproc ({seqproc_name}): {len(seqproc_ids):>12,} emit IDs")
    print(f"  matchbox:                  {len(matchbox_ids):>12,} emit IDs")
    print(f"  splitcode:                 {len(splitcode_ids):>12,} emit IDs")

    if not (seqproc_ids and matchbox_ids and splitcode_ids):
        print("\n[WARN] At least one tool's emit-ID file is missing under")
        print(f"       {CONCORDANCE_DIR}/")
        print("       Run scripts/concordance_analysis.py --datasets splitseq_pe")
        print("       first to populate the ID caches, then re-run this script.")
        sys.exit(2)

    # ---- Step 3: compute intersections ----
    print("\n[3/4] Computing intersections with V_total...")
    tools = [
        ("seqproc",    seqproc_ids),
        ("matchbox",   matchbox_ids),
        ("splitcode",  splitcode_ids),
    ]
    out = {
        "v_total": len(v_total),
        "total_input_reads": 86_820_578,
        "validator": "SplitSeqValidityAnalyzer "
                     "(paired-end, str.find linker on R2, Hamming<=1 whitelist)",
        "tools": {},
    }
    print()
    print(f"{'Tool':<14} {'Emit':>14} {'∩ V_total':>14} {'Recall %':>10}  "
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
        print(f"{label:<14} {len(ids):>14,} {len(inter):>14,} {recall:>9.2f}%  "
              f"{len(missing):>20,}")

    # ---- Step 4: save JSON ----
    out_json = outdir / "pe_recall_vtotal_results.json"
    with open(out_json, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {out_json}")
    print("\nThese are the EXACT numerator/denominator pairs for the SPLiT-seq PE")
    print("row of Table 3 under the V_total-as-denominator framing:")
    print(f"  Denominator (V_total): {len(v_total):,}")
    print(f"  Numerator (tool ∩ V_total): see table above per tool.")


if __name__ == "__main__":
    main()
