#!/usr/bin/env python3
"""
Synthesize a 10x Chromium Flex-like multiplexed dataset.

R1 layout (28 bp): [16 bp GEM cell barcode][12 bp UMI]
R2 layout (86 bp): [50 bp cDNA insert][16 bp constant][12 bp pCS][8 bp probe barcode]

We simulate 4 biological samples with known probe barcodes and fixed mixing
proportions, drawing cell barcodes from a pool of 500 GEMs. The resulting
FASTQ is piped into seqproc with  1{b[16]u[12]}2{r:x[28]s[8]}  to extract
cell, UMI, and sample barcodes. Downstream analysis verifies:

  1. We recover the expected 4 probe barcodes (and only those) as unique s intervals.
  2. Observed probe-barcode proportions match the simulation within sampling error.
  3. Cell-barcode x sample-barcode cross-tab is dense and clean (each cell only
     appears with a single sample, consistent with a true multiplex design).
"""
import random
import sys

random.seed(42)

N_READS = 50_000
N_CELLS = 500
SAMPLES = {
    "SAMPLE01": "ACGTACGT",
    "SAMPLE02": "TGCATGCA",
    "SAMPLE03": "GATTACAG",
    "SAMPLE04": "CCGGCCGG",
}
SAMPLE_NAMES = list(SAMPLES.keys())
SAMPLE_FRACS = [0.40, 0.30, 0.20, 0.10]
CONSTANT = "CAGAGCAATACGACTCACTATAGGG"[:16]
PCS = "TTTCTTATATGGG"[:12]

def rand_seq(n):
    return "".join(random.choices("ACGT", k=n))

# pre-generate cell barcodes, UMIs per cell get generated on the fly
cells = [rand_seq(16) for _ in range(N_CELLS)]

# assign each cell to a unique sample so cells and samples align cleanly
cell_to_sample_idx = [
    random.choices(range(len(SAMPLES)), weights=SAMPLE_FRACS, k=1)[0]
    for _ in range(N_CELLS)
]

r1_path = sys.argv[1]
r2_path = sys.argv[2]

with open(r1_path, "w") as r1_fh, open(r2_path, "w") as r2_fh:
    for i in range(N_READS):
        cell_idx = random.randrange(N_CELLS)
        cell_bc = cells[cell_idx]
        umi = rand_seq(12)
        sample_idx = cell_to_sample_idx[cell_idx]
        sample_bc = SAMPLES[SAMPLE_NAMES[sample_idx]]

        cdna = rand_seq(50)
        r1_seq = cell_bc + umi
        r2_seq = cdna + CONSTANT + PCS + sample_bc
        assert len(r1_seq) == 28
        assert len(r2_seq) == 86

        q1 = "I" * len(r1_seq)
        q2 = "I" * len(r2_seq)
        read_id = f"synth.flex.{i}"

        r1_fh.write(f"@{read_id} 1/1\n{r1_seq}\n+\n{q1}\n")
        r2_fh.write(f"@{read_id} 1/2\n{r2_seq}\n+\n{q2}\n")

# ground truth written next to data
with open(sys.argv[3], "w") as truth:
    truth.write("# cell_barcode\tsample_name\tsample_barcode\n")
    for cb, s_idx in zip(cells, cell_to_sample_idx):
        truth.write(f"{cb}\t{SAMPLE_NAMES[s_idx]}\t{SAMPLES[SAMPLE_NAMES[s_idx]]}\n")

print(f"wrote {N_READS} reads across {N_CELLS} cells, 4 samples")
print(f"expected sample fractions: {dict(zip(SAMPLE_NAMES, SAMPLE_FRACS))}")
