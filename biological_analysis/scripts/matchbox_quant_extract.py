#!/usr/bin/env python3
"""Stitch matchbox TSV (id, umi, bc3, bc2, rest) into a symmetric (cDNA, 34bp barcode) pair.
barcode read = umi(10) + bc3(8) + bc2(8) + rest[:8]; observed, matching seqproc/splitcode. Usage:

  matchbox -s biological_analysis/configs/splitseq_matchbox.mb -t 8 R2.fq \
    | python3 matchbox_quant_extract.py --r1 R1.fq --out-cdna mb_cdna.fq --out-bc mb_bc.fq
"""
import argparse, sys

def read_fastq(fh):
    while True:
        h = fh.readline()
        if not h: break
        s = fh.readline(); fh.readline(); q = fh.readline()
        yield h.split()[0][1:], s.rstrip("\n"), q.rstrip("\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1", required=True)
    ap.add_argument("--out-cdna", required=True)
    ap.add_argument("--out-bc", required=True)
    args = ap.parse_args()
    bc = {}
    for line in sys.stdin:
        p = line.rstrip("\n").split("\t")
        if len(p) != 5: continue
        rid, umi, b3, b2, rest = p
        if len(umi) == 10 and len(b3) == 8 and len(b2) == 8 and len(rest) >= 8:
            bc[rid] = umi + b3 + b2 + rest[:8]          # 34bp
    n = 0
    with open(args.r1) as r1, open(args.out_cdna, "w") as oc, open(args.out_bc, "w") as ob:
        for rid, s, q in read_fastq(r1):
            if rid in bc:
                oc.write(f"@{rid}\n{s}\n+\n{q}\n")
                ob.write(f"@{rid}\n{bc[rid]}\n+\n{'I'*34}\n")
                n += 1
    sys.stderr.write(f"wrote {n} paired (cDNA, 34bp barcode) reads\n")

if __name__ == "__main__":
    main()
