#!/usr/bin/env python3
"""External-truth precision/recall for simulated LR reads, WITH barcode correction.

Truth is in each read id: sim_{i}_bc1_{b1}_bc2_{b2}_bc3_{b3}. Tools emit the RAW
extracted barcode (still carrying read errors), so we snap each to its nearest
whitelist entry within Hamming 1 (trying the reverse complement too, for reads the
tool reported in the opposite orientation) before comparing, mirroring the whitelist
correction the tools/downstream apply. A read is CORRECT iff all three corrected
barcodes equal the truth.

  precision = correct / emitted        recall = correct / n_simulated

Input: TSV `read_id  bc1  bc2  bc3` (matchbox emits this; adapters convert seqproc /
splitcode output to the same columns).

  score_sim_precision.py <assignments.tsv> <n_simulated>
"""
import sys, os, re

TRUTH = re.compile(r"sim_\d+_bc1_([ACGT]+)_bc2_([ACGT]+)_bc3_([ACGT]+)")
_RC = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(_RC)[::-1]


def load(p):
    return [l.strip() for l in open(p) if l.strip()]


def snap_dict(entries):
    """map exact + every Hamming-1 variant -> its canonical whitelist entry."""
    m = {}
    for e in entries:
        m.setdefault(e, e)
        for i in range(len(e)):
            for c in "ACGT":
                m.setdefault(e[:i] + c + e[i + 1:], e)
    return m


def main():
    tsv, n_sim = sys.argv[1], int(sys.argv[2])
    wl = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "configs", "seqproc")
    m23 = snap_dict(load(os.path.join(wl, "splitseq_bc23_whitelist.txt")))
    m1 = snap_dict(load(os.path.join(wl, "splitseq_bc1_whitelist_6bp.txt")))

    def snap(b, m):
        return m.get(b) or m.get(revcomp(b))

    emitted = correct = malformed = 0
    for line in open(tsv):
        p = line.rstrip("\n").split("\t")
        if len(p) < 4:
            malformed += 1
            continue
        rid, rb1, rb2, rb3 = p[0], p[1], p[2], p[3]
        t = TRUTH.match(rid)
        if not t:
            malformed += 1
            continue
        emitted += 1
        c1, c2, c3 = snap(rb1, m1), snap(rb2, m23), snap(rb3, m23)
        if (c1, c2, c3) == (t.group(1), t.group(2), t.group(3)):
            correct += 1
    P = 100 * correct / emitted if emitted else 0.0
    R = 100 * correct / n_sim if n_sim else 0.0
    F1 = 2 * P * R / (P + R) if (P + R) else 0.0
    print(f"emitted={emitted}  correct={correct}  malformed={malformed}")
    print(f"Precision={P:.3f}  Recall={R:.3f}  F1={F1 / 100:.3f}")


if __name__ == "__main__":
    main()
