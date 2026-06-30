#!/usr/bin/env python3
"""Build symmetric (cDNA, barcode-read) FASTQs from splitcode -x extraction.

splitcode emits 4 consecutive --x-only records per read (umi10, bc3-8, bc2-8, bc1-8).
We stitch them into a 34bp observed barcode read (UMI 10 + CB 24) and pair it with the
original cDNA (R1) for reads with a complete barcode. Output matches the seqproc observed
geom byte-for-byte; STARsolo then does all barcode correction. Usage:

  splitcode -c biological_analysis/configs/splitseq_extract.config -N 2 -t 8 --x-only -p \
    -x '1:0<u[10]>,<b3[8]>{linker1},{linker1}<b2[8]>,{linker2}<b1[8]>' R1.fq R2.fq \
    | python3 splitcode_quant_extract.py --r1 R1.fq --out-cdna sc_cdna.fq --out-bc sc_bc.fq
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
    ap.add_argument("--r1", required=True, help="original cDNA reads (R1)")
    ap.add_argument("--out-cdna", required=True)
    ap.add_argument("--out-bc", required=True)
    ap.add_argument("--xstream", default="-", help="splitcode --x-only stream (default stdin)")
    args = ap.parse_args()

    # collect 4 extracted segments per read id from the x-only stream
    seg = {}
    order = []
    xs = sys.stdin if args.xstream == "-" else open(args.xstream)
    for rid, s, q in read_fastq(xs):
        if rid not in seg:
            seg[rid] = []; order.append(rid)
        seg[rid].append((s, q))

    # stitch complete reads: umi(10)+bc3(8)+bc2(8)+bc1(8)=34
    bc = {}
    for rid in order:
        segs = seg[rid]
        if len(segs) == 4 and len(segs[0][0]) == 10 and all(len(x[0]) == 8 for x in segs[1:]):
            bc[rid] = ("".join(x[0] for x in segs), "".join(x[1] for x in segs))

    # pair with cDNA, preserving R1 order
    n = 0
    with open(args.r1) as r1, open(args.out_cdna, "w") as oc, open(args.out_bc, "w") as ob:
        for rid, s, q in read_fastq(r1):
            if rid in bc:
                bseq, bq = bc[rid]
                oc.write(f"@{rid}\n{s}\n+\n{q}\n")
                ob.write(f"@{rid}\n{bseq}\n+\n{bq}\n")
                n += 1
    sys.stderr.write(f"wrote {n} paired (cDNA, 34bp barcode) reads\n")

if __name__ == "__main__":
    main()
