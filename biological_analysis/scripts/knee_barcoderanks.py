#!/usr/bin/env python3
"""Barcode-rank knee/inflection via a faithful Python port of DropletUtils::barcodeRanks
(Lun et al. 2019, EmptyDrops). Used because the R package's HDF5 stack (Rhdf5lib -> libsz) was not
available on the cluster; the algorithm is the reference one.

barcodeRanks fits a smoothing spline to the log-log curve of total UMIs vs barcode rank (over
barcodes with total > lower), then reports two points:
  knee       = maximum curvature (the corner at the top of the plateau),
  inflection = steepest descent, i.e. the top of the cell/empty cliff.
The inflection is the feature that coincides with the called-cell boundary and is what we report; on
synthetic curves with a known cliff it is exact (see test_knee_point.py). We use scipy's GCV
smoothing spline in place of R's smooth.spline; the inflection is insensitive to that choice.

  knee_barcoderanks.py <name1>:<Gene_dir1> <name2>:<Gene_dir2> ...
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from scipy.interpolate import make_smoothing_spline
from concordance_helpers import load_star_raw, per_barcode_umi

def barcode_ranks(totals, lower=100):
    """Port of DropletUtils::barcodeRanks. Returns knee and inflection as UMI thresholds and ranks."""
    totals = np.asarray(totals, float); totals = totals[totals > 0]
    o = np.sort(totals)[::-1]
    change = np.concatenate([[True], o[1:] != o[:-1]])            # run-length encode ties
    idx = np.flatnonzero(change)
    counts = np.diff(np.concatenate([idx, [len(o)]]))
    vals = o[idx]
    run_rank = np.cumsum(counts) - (counts - 1) / 2.0             # averaged rank per run
    keep = vals > lower
    if keep.sum() < 5:
        return None
    x = np.log10(run_rank[keep]); y = np.log10(vals[keep])
    order = np.argsort(x); x, y = x[order], y[order]
    spl = make_smoothing_spline(x, y)                            # GCV smoothing spline (~smooth.spline)
    d1 = spl.derivative(1)(x); d2 = spl.derivative(2)(x)
    curv = d2 / (1 + d1 ** 2) ** 1.5
    knee_umi = 10 ** y[int(np.argmin(curv))]                     # most negative curvature
    infl_umi = 10 ** y[int(np.argmin(d1))]                       # steepest descent (the cliff)
    return {"knee_umi": float(knee_umi), "knee_rank": int((totals >= knee_umi).sum()),
            "infl_umi": float(infl_umi), "infl_rank": int((totals >= infl_umi).sum()),
            "n_barcodes": int(len(totals))}

def main():
    tools = [s.split(":", 1) for s in sys.argv[1:]]
    if not tools:
        sys.exit("usage: knee_barcoderanks.py name:Gene_dir [name:Gene_dir ...]")
    for name, d in tools:
        M, _, _ = load_star_raw(d)
        r = barcode_ranks(per_barcode_umi(M))
        if r is None:
            print(f"{name:<10} (too few barcodes above threshold)"); continue
        print(f"{name:<10} inflection: rank {r['infl_rank']} (UMI>={r['infl_umi']:.0f})   |   "
              f"knee(max-curv): rank {r['knee_rank']} (UMI>={r['knee_umi']:.0f})   |   "
              f"n_barcodes {r['n_barcodes']}")

if __name__ == "__main__":
    main()
