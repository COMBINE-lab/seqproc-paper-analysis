#!/usr/bin/env python3
"""External-truth precision/recall for simulated LR reads.

Truth is encoded in each read id: sim_{i}_bc1_{b1}_bc2_{b2}_bc3_{b3}.
Input is a tool's assignment as TSV `read_id  bc1  bc2  bc3` (matchbox emits this
directly; for seqproc/splitcode, dump read_id + extracted barcodes into the same
3-column form). A read is CORRECT iff all three assigned barcodes equal the truth.

  precision = correct / emitted        (of what the tool output, how much is right)
  recall    = correct / n_simulated    (of all true reads, how many recovered right)

  score_sim_precision.py <tool_assignments.tsv> <n_simulated>
"""
import sys, re

TRUTH = re.compile(r"sim_\d+_bc1_([ACGT]+)_bc2_([ACGT]+)_bc3_([ACGT]+)")


def main():
    tsv, n_sim = sys.argv[1], int(sys.argv[2])
    emitted = correct = malformed = 0
    for line in open(tsv):
        p = line.rstrip("\n").split("\t")
        if len(p) < 4:
            malformed += 1
            continue
        rid, b1, b2, b3 = p[0], p[1], p[2], p[3]
        m = TRUTH.match(rid)
        if not m:
            malformed += 1
            continue
        emitted += 1
        if (b1, b2, b3) == (m.group(1), m.group(2), m.group(3)):
            correct += 1
    P = 100 * correct / emitted if emitted else 0.0
    R = 100 * correct / n_sim if n_sim else 0.0
    F1 = 2 * P * R / (P + R) if (P + R) else 0.0
    print(f"emitted={emitted}  correct={correct}  malformed={malformed}")
    print(f"Precision={P:.3f}  Recall={R:.3f}  F1={F1/100:.3f}")


if __name__ == "__main__":
    main()
