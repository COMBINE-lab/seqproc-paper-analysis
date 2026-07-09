#!/usr/bin/env python3
"""Long-read precision against the short-read real-cell set (external truth).

real_cells.txt (from short_read_cb_caller.py, built on the clean short reads) is the set of
true cells. A long-read tool's emitted read is CORRECT iff its (whitelist-snapped) cell
barcode bc1_bc2_bc3 is in that set. Uncorrectable / non-real barcodes count as incorrect.

  precision = correct / emitted

Input formats:
  matchbox : read_id\tbc1\tbc2\tbc3\t...   (splitseq_singleend_dual.mb emits this)
  seqproc  : fastq whose seq is umi(10) bc3(8) bc2(8) bc1(8); we take the 3 barcodes
             (bc1 first 6 of its 8) from fixed offsets

  score_lr_precision.py real_cells.txt tool_output --format {matchbox,seqproc}
"""
import sys, os, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from short_read_cb_caller import snap_dict, load


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("real_cells")
    ap.add_argument("tool_output")
    ap.add_argument("--format", choices=["matchbox", "seqproc"], default="matchbox")
    a = ap.parse_args()

    wl = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "configs", "seqproc")
    m23 = snap_dict(load(os.path.join(wl, "splitseq_bc23_whitelist.txt")))
    m1 = snap_dict(load(os.path.join(wl, "splitseq_bc1_whitelist_6bp.txt")))
    real = set(l.strip() for l in open(a.real_cells) if l.strip())

    def records():
        if a.format == "matchbox":
            for line in open(a.tool_output):
                p = line.rstrip("\n").split("\t")
                if len(p) >= 4:
                    yield p[1], p[2], p[3]
        else:  # seqproc fastq: umi[0:10] bc3[10:18] bc2[18:26] bc1[26:34] (bc1 first 6)
            with open(a.tool_output) as f:
                while True:
                    h = f.readline()
                    if not h:
                        break
                    s = f.readline().strip(); f.readline(); f.readline()
                    if len(s) >= 32:
                        yield s[26:32], s[18:26], s[10:18]   # bc1(6), bc2, bc3

    emitted = correct = 0
    for rb1, rb2, rb3 in records():
        emitted += 1
        b1, b2, b3 = m1.get(rb1), m23.get(rb2), m23.get(rb3)
        if b1 and b2 and b3 and f"{b1}_{b2}_{b3}" in real:
            correct += 1
    P = 100 * correct / emitted if emitted else 0.0
    print(f"emitted={emitted}  correct={correct}  precision={P:.3f}%")


if __name__ == "__main__":
    main()
