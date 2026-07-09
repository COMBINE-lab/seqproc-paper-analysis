#!/usr/bin/env python3
"""Short-read cell-barcode caller -> the external-truth REAL-CELL set for LR precision.

From the clean barcode read (R2), locate linker1 and linker2, read bc3/bc2/bc1 at
linker-relative positions, snap each to its whitelist (Hamming 1), form the cell barcode
CB = bc1_bc2_bc3, count CBs across reads, and emit the set of REAL cells (CB count >=
--min-count). Because short reads are low-error, this set is a reliable, tool-independent
reference: a long-read tool's assignment is "correct" iff its CB is in this set.

  short_read_cb_caller.py R2.fastq --chem lr --out real_cells.txt [--min-count 10] [--sample N]
"""
import sys, os, argparse
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import edit_tolerant_validity as V


def load(p):
    return [l.strip() for l in open(p) if l.strip()]


def snap_dict(entries):
    m = {}
    for e in entries:
        m.setdefault(e, e)
        for i in range(len(e)):
            for c in "ACGT":
                m.setdefault(e[:i] + c + e[i + 1:], e)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fastq")
    ap.add_argument("--chem", choices=["pe", "lr"], default="lr")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-count", type=int, default=10)
    ap.add_argument("--sample", type=int, default=0)
    ap.add_argument("--max-linker-edit", type=int, default=3)   # short reads are clean
    a = ap.parse_args()

    L1, L2 = V.LINKERS[a.chem]
    wl = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "configs", "seqproc")
    wl23 = load(os.path.join(wl, "splitseq_bc23_whitelist.txt")); n23 = len(wl23[0])
    wl1 = load(os.path.join(wl, "splitseq_bc1_whitelist_6bp.txt")); n1 = len(wl1[0])
    m23, m1 = snap_dict(wl23), snap_dict(wl1)

    import gzip
    opener = gzip.open if a.fastq.endswith(".gz") else open
    cb = Counter(); total = seen = 0
    with opener(a.fastq, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                break
            seq = f.readline().strip(); f.readline(); f.readline()
            if a.sample and total >= a.sample:
                break
            total += 1
            l1 = V.find(L1, seq, a.max_linker_edit)
            if not l1:
                continue
            s1, e1 = l1
            bc3 = m23.get(seq[s1 - n23:s1])
            bc2 = m23.get(seq[e1 + 1:e1 + 1 + n23])
            if not bc3 or not bc2:
                continue
            l2 = V.find(L2, seq[e1 + 1:], a.max_linker_edit)
            if not l2:
                continue
            e2 = l2[1] + e1 + 1
            bc1 = m1.get(seq[e2 + 1:e2 + 1 + n1])
            if not bc1:
                continue
            cb[(bc1, bc2, bc3)] += 1
            seen += 1

    real = {c for c, n in cb.items() if n >= a.min_count}
    with open(a.out, "w") as o:
        for b1, b2, b3 in sorted(real):
            o.write(f"{b1}_{b2}_{b3}\n")
    print(f"reads={total}  with_valid_CB={seen}  distinct_CBs={len(cb)}  "
          f"real_cells(count>={a.min_count})={len(real)}  wrote {a.out}")


if __name__ == "__main__":
    main()
