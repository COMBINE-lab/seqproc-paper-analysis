#!/usr/bin/env python3
"""Simulate LR-SPLiT-seq reads with KNOWN barcodes + PacBio-like error.

External-truth precision benchmark: each read header encodes its true bc1/bc2/bc3,
so after running a tool, precision = (emitted reads whose assigned barcode == truth)
/ (emitted reads). The truth is the simulated barcode, independent of every tool, so
this is non-circular (unlike GT_edit, which uses edit-distance matching like the tools).

Read layout (forward): [5' junk][UMI 10][bc3 8][L1 30][bc2 8][L2 22][bc1 6][cDNA];
half the reads are reverse-complemented (PacBio random orientation). Substitution /
insertion / deletion errors are applied per base to mimic the degraded long-read library.

  simulate_lr_reads.py --n 100000 --out sim.fastq [--sub 0.04 --ins 0.03 --del 0.03]
"""
import os, random, argparse

L1 = "GTGGCCGATGTTTCGCATCGGCGTACGACT"   # LR linker 1 (30 bp)
L2 = "ATCCACGTGCTTGAGACTGTGG"           # LR linker 2 (22 bp)
_RC = str.maketrans("ACGT", "TGCA")


def load(p):
    return [l.strip() for l in open(p) if l.strip()]


def rand_seq(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def add_errors(seq, sub, ins, dele):
    out = []
    for b in seq:
        r = random.random()
        if r < dele:
            continue
        if r < dele + ins:
            out.append(random.choice("ACGT"))
            out.append(b)
        elif r < dele + ins + sub:
            out.append(random.choice([x for x in "ACGT" if x != b]))
        else:
            out.append(b)
    return "".join(out)


def revcomp(s):
    return s.translate(_RC)[::-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=100000)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sub", type=float, default=0.04)
    ap.add_argument("--ins", type=float, default=0.03)
    ap.add_argument("--del", dest="dele", type=float, default=0.03)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    random.seed(a.seed)

    wl = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "configs", "seqproc")
    bc23 = load(os.path.join(wl, "splitseq_bc23_whitelist.txt"))
    bc1s = load(os.path.join(wl, "splitseq_bc1_whitelist_6bp.txt"))

    with open(a.out, "w") as o:
        for i in range(a.n):
            b1, b2, b3 = random.choice(bc1s), random.choice(bc23), random.choice(bc23)
            umi = rand_seq(10)
            core = (rand_seq(random.randint(0, 20)) + umi + b3 + L1 + b2 + L2 + b1
                    + rand_seq(random.randint(100, 300)))
            read = add_errors(core, a.sub, a.ins, a.dele)
            if random.random() < 0.5:
                read = revcomp(read)
            rid = f"sim_{i}_bc1_{b1}_bc2_{b2}_bc3_{b3}"        # truth in the header
            o.write(f"@{rid}\n{read}\n+\n{'I' * len(read)}\n")
    print(f"wrote {a.n} reads to {a.out}  (sub={a.sub} ins={a.ins} del={a.dele})")


if __name__ == "__main__":
    main()
