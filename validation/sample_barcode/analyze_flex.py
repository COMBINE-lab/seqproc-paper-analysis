#!/usr/bin/env python3
"""
Downstream validation of seqproc `s` (sample barcode) feature on synthetic
10x Chromium Flex multiplexed data.

We extract cell barcode from R1[:16] and sample barcode from the trailing 8 bp
of each processed R2 read, then cross-check against the ground-truth
cell -> sample_name mapping.

Checks:
  1. We recover exactly the 4 known probe barcodes.
  2. Observed sample proportions match the simulated 40/30/20/10 mix within
     sampling tolerance.
  3. Every cell is paired with its true sample in every read (no mixing).
"""
import argparse
import collections
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_WORK = os.environ.get("FLEX_WORK_DIR", os.path.join(_HERE, "work"))

_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--work-dir", default=_DEFAULT_WORK,
                 help="directory with flex_truth.tsv and out_flex_R{1,2}.fastq "
                      "(default: ./work or $FLEX_WORK_DIR)")
_args, _ = _ap.parse_known_args()

TRUTH_TSV = os.path.join(_args.work_dir, "flex_truth.tsv")
R1_FQ = os.path.join(_args.work_dir, "out_flex_R1.fastq")
R2_FQ = os.path.join(_args.work_dir, "out_flex_R2.fastq")

EXPECTED = {
    "ACGTACGT": ("SAMPLE01", 0.40),
    "TGCATGCA": ("SAMPLE02", 0.30),
    "GATTACAG": ("SAMPLE03", 0.20),
    "CCGGCCGG": ("SAMPLE04", 0.10),
}

def parse_truth(path):
    mapping = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cb, name, sbc = line.rstrip("\n").split("\t")
            mapping[cb] = (name, sbc)
    return mapping

def iter_seqs(path):
    with open(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline().rstrip("\n")
            _ = fh.readline()
            _ = fh.readline()
            yield seq

def main():
    truth = parse_truth(TRUTH_TSV)
    print(f"loaded truth for {len(truth)} cells")

    # expected read fractions are the cell-assignment fractions, because each
    # read picks a cell uniformly. The nominal 40/30/20/10 weights are the
    # *prior* on cell assignment, not the realized cell distribution.
    truth_cell_counts = collections.Counter(name for name, _ in truth.values())
    realized_fracs = {
        sbc: truth_cell_counts[name] / len(truth)
        for sbc, (name, _) in EXPECTED.items()
    }

    sample_counts = collections.Counter()
    cross_tab = collections.defaultdict(lambda: collections.Counter())
    mismatches = 0
    total = 0

    for r1, r2 in zip(iter_seqs(R1_FQ), iter_seqs(R2_FQ)):
        cell_bc = r1[:16]
        sample_bc = r2[-8:]
        sample_counts[sample_bc] += 1
        cross_tab[cell_bc][sample_bc] += 1
        true_name, true_sbc = truth.get(cell_bc, (None, None))
        if sample_bc != true_sbc:
            mismatches += 1
        total += 1

    print(f"\n== check 1: unique probe barcodes ==")
    recovered = set(sample_counts)
    expected_set = set(EXPECTED)
    print(f"  recovered = {sorted(recovered)}")
    print(f"  expected  = {sorted(expected_set)}")
    print(f"  match     = {recovered == expected_set}")
    assert recovered == expected_set, "probe-barcode set mismatch"

    print(f"\n== check 2: sample proportions ==")
    print(f"  ({'nominal':7s} = prior on cell assignment; "
          f"{'realized':8s} = fraction of the 500 cells actually drawn)")
    print(f"  {'sample':10s} {'barcode':10s} {'nominal':>8s} "
          f"{'realized':>9s} {'observed':>9s} {'|obs-real|':>10s}")
    for sbc, (name, nominal) in EXPECTED.items():
        obs = sample_counts[sbc] / total
        real = realized_fracs[sbc]
        diff = abs(obs - real)
        tag = "OK" if diff < 0.01 else "WARN"
        print(f"  {name:10s} {sbc:10s} {nominal:>8.4f} "
              f"{real:>9.4f} {obs:>9.4f} {diff:>10.4f}  [{tag}]")
        assert diff < 0.01, f"{name} observed vs realized off by {diff:.4f}"

    print(f"\n== check 3: cell x sample consistency ==")
    print(f"  total reads scanned          = {total}")
    print(f"  reads w/ sample != truth(cb) = {mismatches}")
    assert mismatches == 0, f"{mismatches} reads disagree with truth"

    n_cells = len(cross_tab)
    single_sample_cells = sum(1 for c in cross_tab.values() if len(c) == 1)
    multi_sample_cells = n_cells - single_sample_cells
    print(f"  cells observed               = {n_cells}")
    print(f"  cells with single sample     = {single_sample_cells}")
    print(f"  cells with >1 sample (leak)  = {multi_sample_cells}")
    assert multi_sample_cells == 0, "cell barcode appeared with multiple samples"

    print("\n[PASS] synthetic Flex recovery is clean: 4/4 probe barcodes, "
          "proportions within 1%, no cell-sample leak.")

if __name__ == "__main__":
    main()
