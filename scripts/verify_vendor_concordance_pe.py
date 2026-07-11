#!/usr/bin/env python3
"""
Verify the SPLiT-seq PE vendor-concordance numbers reported in the paper,
WITHOUT needing split-pipe installed.

split-pipe is proprietary (Parse Biosciences); its valid-barcode read-ID set
on the 10M subset (SRR6750041) is shipped as
    results_final/splitpipe_valid_ids_10M.txt.gz
(just read IDs, no sequence). This script intersects that ground-truth set
with seqproc's and splitcode's emitted reads and reports precision / recall /
F1 / Jaccard, then checks them against the committed expected values.

You provide the two tools' R2 outputs, produced by running them on the SAME
10M subset with the committed configs:

  seqproc:
    seqproc --geom configs/seqproc/splitseq_filter_edit.geom \\
      --file1 R1_10M.fastq --file2 R2_10M.fastq \\
      --out1 sp_R1.fq --out2 sp_R2.fq --threads 8 \\
      -a configs/seqproc/splitseq_bc3_seq2seq.tsv \\
      -a configs/seqproc/splitseq_bc2_seq2seq.tsv \\
      -a configs/seqproc/splitseq_bc1_seq2seq.tsv

  splitcode:
    splitcode -c configs/splitcode/splitseq_paper.config --assign -N 2 -t 8 \\
      -m mapping.txt -o sc_R1.fq,sc_R2.fq R1_10M.fastq R2_10M.fastq

Then:
    python3 scripts/verify_vendor_concordance_pe.py \\
        --seqproc-r2 sp_R2.fq --splitcode-r2 sc_R2.fq
"""
import argparse
import gzip
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
GT_GZ = HERE / "results_final" / "splitpipe_valid_ids_10M.txt.gz"
EXPECTED = HERE / "results_final" / "vendor_concordance_pe_10M.json"


def fastq_ids(path):
    s = set()
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        for i, l in enumerate(f):
            if i % 4 == 0:
                s.add(l.split()[0].lstrip("@"))
    return s


def load_gt():
    s = set()
    with gzip.open(GT_GZ, "rt") as f:
        for l in f:
            l = l.strip()
            if l:
                s.add(l)
    return s


def metrics(tool, gt):
    inter = len(tool & gt)
    P = inter / len(tool) if tool else 0
    R = inter / len(gt) if gt else 0
    F1 = 2 * P * R / (P + R) if (P + R) else 0
    J = inter / len(tool | gt) if (tool or gt) else 0
    return inter, 100 * P, 100 * R, F1, J


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seqproc-r2", required=True, help="seqproc R2 output on the 10M subset")
    ap.add_argument("--splitcode-r2", required=True, help="splitcode R2 output on the 10M subset")
    args = ap.parse_args()

    if not GT_GZ.exists():
        raise SystemExit(f"missing ground-truth IDs: {GT_GZ}")

    print(f"Loading split-pipe ground-truth IDs from {GT_GZ.name} ...")
    gt = load_gt()
    print(f"  split-pipe valid: {len(gt):,}\n")

    exp = json.load(open(EXPECTED)) if EXPECTED.exists() else {"tools": {}}
    rows = [("seqproc", args.seqproc_r2), ("splitcode", args.splitcode_r2)]

    print(f"{'Tool':<11}{'Emit':>12}{'Inter':>12}{'Prec%':>8}{'Recall%':>9}{'F1':>8}{'Jacc':>7}  check")
    print("-" * 76)
    ok = True
    for tool, path in rows:
        ids = fastq_ids(path)
        inter, P, R, F1, J = metrics(ids, gt)
        e = exp.get("tools", {}).get(tool, {})
        # tolerance: 0.2 percentage points on precision/recall
        check = "n/a"
        if e:
            dp = abs(P - e.get("precision", P))
            dr = abs(R - e.get("recall", R))
            passed = dp <= 0.2 and dr <= 0.2
            ok = ok and passed
            check = "PASS" if passed else f"MISMATCH (exp P={e['precision']} R={e['recall']})"
        print(f"{tool:<11}{len(ids):>12,}{inter:>12,}{P:>8.2f}{R:>9.2f}{F1:>8.4f}{J:>7.3f}  {check}")

    print()
    if ok:
        print("OK -- measured numbers match the committed expected values.")
    else:
        print("WARNING -- numbers differ from committed expected values; investigate before publishing.")


if __name__ == "__main__":
    main()
