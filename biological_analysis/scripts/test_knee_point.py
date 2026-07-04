#!/usr/bin/env python3
"""Ground-truth tests for the barcode-rank knee detector (count_concordance.knee_point).

knee_point returns the inflection (steepest descent) of the log-log barcode-rank curve, i.e. the
top of the cell/empty cliff. We validate it against synthetic curves whose true cliff rank is KNOWN
by construction: a mixture of high-UMI cells over low-UMI ambient barcodes, where the cliff sits at
~n_cells. This is exactly the shape of real barcode-rank data.

Scope (important, and asserted honestly): this detector is validated for realistic barcode-rank
curves, which are bimodal with a sharp cell/empty cliff. It is NOT a general-purpose knee detector.
Broad, gradual transitions (e.g. a wide sigmoid whose inflection is at high rank) are out of scope
and can overshoot; real barcode-rank data is not shaped that way.

  python test_knee_point.py     # prints PASS/FAIL per test
  pytest test_knee_point.py
"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from count_concordance import knee_point

def _mixture(n_cells, sigma_c=0.6, amb=3.0, sigma_a=0.9, n=15000, noise=0.0, seed=0):
    """Cells (high UMI) over ambient (low UMI); the true cliff is at ~n_cells."""
    rng = np.random.default_rng(seed)
    umi = np.concatenate([rng.lognormal(np.log(4000), sigma_c, n_cells),
                          rng.lognormal(np.log(amb), sigma_a, n - n_cells)])
    if noise:
        umi = umi * np.exp(rng.normal(0, noise, n))
    return np.sort(np.clip(np.round(umi), 1, None))[::-1]

def _knee(umi, **k):
    return knee_point(np.asarray(umi, float), **k)[0]

def test_accuracy_on_known_cliff():
    """Detected knee lands within 15% of the true cliff across cell count and cliff sharpness."""
    for n_cells in (100, 220, 500, 1000, 2000):
        for sigma_c in (0.4, 0.6, 0.8, 1.0):
            ratio = np.median([_knee(_mixture(n_cells, sigma_c=sigma_c, seed=s)) / n_cells
                               for s in range(5)])
            assert 0.85 <= ratio <= 1.15, f"n_cells={n_cells} sigma={sigma_c}: ratio {ratio:.2f}"

def test_floor_stable():
    """umi_floor 5..20 agree (the default 10 sits inside this stable band)."""
    umi = np.asarray(_mixture(220, noise=0.05, seed=1), float)
    vals = {knee_point(umi, umi_floor=f)[0] for f in (5, 10, 20)}
    assert max(vals) - min(vals) <= 3, vals

def test_three_tool_consistency():
    """Three independent draws of one distribution agree within ~15 barcodes (like the 3 tools)."""
    for rep in range(8):
        ks = [_knee(_mixture(220, noise=0.05, seed=100 * rep + t)) for t in range(3)]
        assert max(ks) - min(ks) <= 15, ks

def test_edge_cases():
    assert knee_point(np.full(1000, 5.0))[0] is not None            # flat curve, no crash
    assert knee_point(np.array([9., 8, 7, 3, 2, 1])) == (None, None)  # too few barcodes

if __name__ == "__main__":
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn(); print("PASS", name)
            except AssertionError as e:
                fails += 1; print("FAIL", name, "->", e)
    print("ALL PASS" if not fails else f"{fails} FAILED")
    sys.exit(1 if fails else 0)
