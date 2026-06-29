#!/usr/bin/env python3
"""
Keep only read pairs whose barcode read is the canonical fixed length, so the
output is STARsolo-ready (fixed CB and UMI positions). Reads with a barcode
read of any other length are dropped, since a non-canonical length means the
barcode was not extracted at the expected positions and would not form a valid
cell. This is applied to seqproc's output to match splitcode's reformatting
mode, which already emits only fixed-length whitelist-valid barcodes, keeping
the two tools symmetric.

Usage:
  filter_fixed_length_bc.py --cdna R1.fq --barcode R2.fq --out-prefix OUT [--length N]
If --length is omitted, the modal barcode-read length is used.
Writes OUT_cdna.fq and OUT_barcode.fq.
"""
import argparse
from collections import Counter


def read_fastq(path):
    with open(path) as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline()
            p = f.readline()
            q = f.readline()
            yield h, s, p, q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cdna", required=True)
    ap.add_argument("--barcode", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--length", type=int, default=None)
    args = ap.parse_args()

    # determine target length (modal) if not given
    target = args.length
    if target is None:
        c = Counter()
        for _, s, _, _ in read_fastq(args.barcode):
            c[len(s.strip())] += 1
        target = c.most_common(1)[0][0]
        print(f"  modal barcode length = {target}  (distribution: {dict(c.most_common())})")

    kept = 0
    total = 0
    with open(f"{args.out_prefix}_cdna.fq", "w") as oc, open(f"{args.out_prefix}_barcode.fq", "w") as ob:
        for (h1, s1, p1, q1), (h2, s2, p2, q2) in zip(read_fastq(args.cdna), read_fastq(args.barcode)):
            total += 1
            if len(s2.strip()) == target:
                oc.write(h1 + s1 + p1 + q1)
                ob.write(h2 + s2 + p2 + q2)
                kept += 1
    print(f"  kept {kept:,} / {total:,} pairs ({100*kept/total:.1f}%) at length {target}")


if __name__ == "__main__":
    main()
