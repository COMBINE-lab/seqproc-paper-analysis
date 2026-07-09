#!/usr/bin/env python3
"""Edit-tolerant genuine-barcode validity for the Valid/precision column.

Builds the set of read IDs that genuinely carry the SPLiT-seq PE barcode
structure while ALLOWING indels, unlike the strict V_total (fixed offset,
Hamming linker). Linker1 and Linker2 are each located by edlib infix search
within edit distance <= --max-linker-edit, and the three barcodes, read at
positions relative to the located linkers, are each required within Hamming
distance 1 of their whitelist. Because it locates the linker instead of
assuming a fixed offset, it credits reads whose linker carries an indel (the
recoveries seqproc/matchbox are built for) while still failing reads that lack
a linker entirely (splitcode's over-emission). It is therefore an unbiased
precision reference.

Precision for a tool = |emitted_ids & valid_ids| / |emitted_ids|.

Usage:
  edit_tolerant_validity.py R2.fastq --out valid_ids.txt [--sample N] [--max-linker-edit 6]
"""
import sys, os, argparse, json
import edlib

# Linker sequences per chemistry. PE and LR-SPLiT-seq differ (LR L1 has an A at
# position 8; LR L2 is the shorter 22 bp variant).
LINKERS = {
    "pe": ("GTGGCCGCTGTTTCGCATCGGCGTACGACT", "ATCCACGTGCTTGAGAGGCCAGAGCATTCG"),
    "lr": ("GTGGCCGATGTTTCGCATCGGCGTACGACT", "ATCCACGTGCTTGAGACTGTGG"),
}
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def ham1_set(path):
    bcs = [l.strip() for l in open(path) if l.strip()]
    s = set(bcs)
    for b in bcs:
        for i in range(len(b)):
            for c in "ACGT":
                if c != b[i]:
                    s.add(b[:i] + c + b[i + 1:])
    return s, (len(bcs[0]) if bcs else 0)


def find(query, target, max_edit):
    r = edlib.align(query, target, mode="HW", task="locations")
    if r["editDistance"] < 0 or r["editDistance"] > max_edit or not r["locations"]:
        return None
    return r["locations"][0]                       # (start, end) inclusive, best match


def genuine(seq, L1, L2, bc23, bc1s, n23, n1, max_edit):
    """True if seq carries linker1 -> bc2 -> linker2 -> bc1 with a whitelist-valid
    bc3 (before L1), bc2 (after L1) and bc1 (after L2), linkers within edit<=max_edit."""
    l1 = find(L1, seq, max_edit)
    if not l1:
        return False
    s1, e1 = l1
    bc3 = seq[s1 - n23:s1]
    bc2 = seq[e1 + 1:e1 + 1 + n23]
    if len(bc3) != n23 or len(bc2) != n23 or bc3 not in bc23 or bc2 not in bc23:
        return False
    l2 = find(L2, seq[e1 + 1:], max_edit)
    if not l2:
        return False
    e2 = l2[1] + e1 + 1
    bc1 = seq[e2 + 1:e2 + 1 + n1]
    return len(bc1) == n1 and bc1 in bc1s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fastq")
    ap.add_argument("--chem", choices=["pe", "lr"], default="pe",
                    help="pe: check forward only. lr: check both orientations (PacBio).")
    ap.add_argument("--out", default=None)
    ap.add_argument("--sample", type=int, default=0)
    ap.add_argument("--max-linker-edit", type=int, default=6)
    a = ap.parse_args()

    L1, L2 = LINKERS[a.chem]
    wl = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "configs", "seqproc")
    bc23, n23 = ham1_set(os.path.join(wl, "splitseq_bc23_whitelist.txt"))
    bc1s, n1 = ham1_set(os.path.join(wl, "splitseq_bc1_whitelist_6bp.txt"))

    valid, total = set(), 0
    with open(a.fastq) as f:
        while True:
            h = f.readline()
            if not h:
                break
            seq = f.readline().strip(); f.readline(); f.readline()
            if a.sample and total >= a.sample:
                break
            total += 1
            rid = h[1:].split()[0]
            ok = genuine(seq, L1, L2, bc23, bc1s, n23, n1, a.max_linker_edit)
            if not ok and a.chem == "lr":
                ok = genuine(revcomp(seq), L1, L2, bc23, bc1s, n23, n1, a.max_linker_edit)
            if ok:
                valid.add(rid)

    print(json.dumps({"fastq": os.path.basename(a.fastq), "chem": a.chem, "total": total,
                      "valid": len(valid),
                      "pct_of_scanned": round(100 * len(valid) / total, 2) if total else 0.0}))
    if a.out:
        with open(a.out, "w") as o:
            o.write("\n".join(sorted(valid)))
        print("wrote", a.out)


if __name__ == "__main__":
    main()
