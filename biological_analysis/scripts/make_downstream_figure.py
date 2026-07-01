#!/usr/bin/env python3
"""Merged downstream-concordance figure (fig:downstream) for the paper, in one 2x2 panel:
(A) barcode-rank knee, (B) per-gene total concordance, (C) joint embedding by tool,
(D) cell-type composition. Reuses the exact logic of count_concordance.py and
biological_analysis.py so the merged figure matches the standalone panels.

  make_downstream_figure.py <outdir> <min_umi> <name1>:<Gene_dir1> <name2>:<dir2> [<name3>:<dir3> ...]
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import matplotlib; matplotlib.use("Agg")
import numpy as np, matplotlib.pyplot as plt
from numpy import log1p, corrcoef
from scipy.stats import spearmanr
from concordance_helpers import load_star_raw, barcode_rank
from biological_analysis import load, process, joint_embedding, MARKERS
from paper_style import set_paper_style, tool_color, panel, save

def main():
    set_paper_style()
    outdir, min_umi = sys.argv[1], int(sys.argv[2])
    tools = [s.split(":", 1) for s in sys.argv[3:]]
    os.makedirs(outdir, exist_ok=True)
    names = [n for n, _ in tools]
    ref = names[0]

    # count-level inputs (knee + per-gene totals)
    M = {}
    for n, d in tools:
        M[n], _, _ = load_star_raw(d)
    pgt = {n: np.asarray(M[n].sum(0)).ravel() for n in names}

    # biological inputs (cell calling, typing, joint embedding)
    proc = {n: process(load(d), min_umi=min_umi) for n, d in tools}
    raw = {n: load(d) for n, d in tools}
    shared = sorted(set.intersection(*[set(proc[n].obs_names) for n in names]))
    fracs = {n: {ct: float((proc[n].obs["cell_type"] == ct).mean()) for ct in MARKERS} for n in names}
    comb, co, _mix = joint_embedding(raw, shared, min_umi)

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    # (A) barcode-rank knee
    ax = axes[0, 0]
    for i, n in enumerate(names):
        r = barcode_rank(M[n]); r = r[r > 0]
        ax.loglog(np.arange(1, len(r) + 1), r, label=f"{n} (n={len(r):,})", color=tool_color(n, i), lw=1.8)
    ax.set_xlabel("Barcode rank"); ax.set_ylabel("Total UMI per barcode")
    ax.set_title("Barcode rank"); ax.legend(loc="lower left")
    ax.grid(True, which="major", alpha=0.25, lw=0.5); panel(ax, "A")

    # (B) per-gene total concordance
    ax = axes[0, 1]; xr = pgt[ref]
    for i, n in enumerate(names[1:], 1):
        y = pgt[n]; m = (xr + y) > 0
        pr = float(corrcoef(log1p(xr[m]), log1p(y[m]))[0, 1]); sp = float(spearmanr(xr[m], y[m]).correlation)
        ax.scatter(log1p(xr[m]), log1p(y[m]), s=6, alpha=0.35, color=tool_color(n, i), edgecolors="none",
                   label=f"{n}:  r={pr:.3f}, ρ={sp:.3f}")
    lim = max(log1p(xr).max(), 1); ax.plot([0, lim], [0, lim], color="0.4", ls="--", lw=1)
    ax.set_xlabel(f"{ref}  log1p(gene UMI)"); ax.set_ylabel("other tool  log1p(gene UMI)")
    ax.set_title("Per-gene total"); ax.legend(loc="upper left", handletextpad=0.2)
    ax.set_aspect("equal", adjustable="box"); panel(ax, "B")

    # (C) joint embedding colored by tool
    ax = axes[1, 0]; u = comb.obsm["X_umap"]
    for i, n in enumerate(names):
        m = comb.obs["tool"].values == n
        ax.scatter(u[m, 0], u[m, 1], s=10, alpha=0.7, label=n, color=tool_color(n, i), edgecolors="none")
    ax.legend(loc="best", markerscale=1.6)
    ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"Joint embedding by tool  (co-clustering = {co:.3f})"); panel(ax, "C")

    # (D) cell-type composition
    ax = axes[1, 1]
    cts = list(MARKERS.keys()); x = np.arange(len(cts)); w = 0.8 / len(names)
    for i, n in enumerate(names):
        ax.bar(x + i * w - 0.4 + w / 2, [fracs[n][c] for c in cts], w, label=n, color=tool_color(n, i))
    ax.set_xticks(x); ax.set_xticklabels(cts, rotation=35, ha="right")
    ax.set_ylabel("Fraction of cells"); ax.set_title("Cell-type composition")
    ax.legend(loc="upper right"); panel(ax, "D")

    fig.tight_layout()
    print("saved", save(fig, os.path.join(outdir, "downstream_concordance")), "(+ .pdf)")

if __name__ == "__main__":
    main()
