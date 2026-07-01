#!/usr/bin/env python3
"""Merged downstream-concordance figure (fig:downstream), one 2x2 panel, NO UMAP/embedding:
(A) barcode-rank knee, (B) per-gene total concordance, (C) quantitative concordance scorecard,
(D) cell-type composition. Panels A/B come from the STARsolo matrices; C/D are read from
biological_metrics.json (written by biological_analysis.py into the same outdir), so run that
first.

Emits TWO versions: `downstream_concordance` (scorecard = agreement, cluster ARI, co-clustering)
and `downstream_concordance_jaccard` (same plus the mean cell-type Jaccard bar).

  make_downstream_figure.py <outdir> <min_umi> <name1>:<Gene_dir1> <name2>:<dir2> [<name3>:<dir3> ...]
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
import matplotlib; matplotlib.use("Agg")
import numpy as np, matplotlib.pyplot as plt
from numpy import log1p, corrcoef
from scipy.stats import spearmanr
from concordance_helpers import load_star_raw, barcode_rank
from paper_style import set_paper_style, tool_color, panel, save

def scorecard(ax, labels, vals, shared_n):
    yb = np.arange(len(labels))[::-1]
    ax.barh(yb, vals, color="#1b7837", height=0.55)
    ax.axvline(1.0, color="0.5", ls="--", lw=1)
    for y_, v in zip(yb, vals):
        if v == v:
            ax.text(min(v, 0.97) - 0.02, y_, f"{v:.3f}", va="center", ha="right", color="white", fontweight="bold")
    ax.set_yticks(yb); ax.set_yticklabels(labels); ax.set_xlim(0, 1.08); ax.set_xlabel("Score (1 = identical)")
    ax.set_title(f"Tool concordance  (n={shared_n} shared cells)")
    ax.spines["left"].set_visible(False)

def main():
    set_paper_style()
    outdir, min_umi = sys.argv[1], int(sys.argv[2])
    tools = [s.split(":", 1) for s in sys.argv[3:]]
    os.makedirs(outdir, exist_ok=True)
    names = [n for n, _ in tools]
    ref = names[0]

    M = {}
    for n, d in tools:
        M[n], _, _ = load_star_raw(d)
    pgt = {n: np.asarray(M[n].sum(0)).ravel() for n in names}

    bm = json.load(open(os.path.join(outdir, "biological_metrics.json"))) if os.path.exists(os.path.join(outdir, "biological_metrics.json")) else {}
    g = lambda k: bm.get(k) if isinstance(bm.get(k), (int, float)) else float("nan")
    shared_n = bm.get("shared_cells", "?")
    fr = bm.get("celltype_fractions", {})

    # scorecard metric sets: default (no mean Jaccard) and the _jaccard variant
    base = [("Cell-type\nagreement", g("celltype_agreement_shared")),
            ("Cluster ARI\n(mean)", g("cluster_ari_mean")),
            ("Co-clustering\n(joint)", g("joint_coclustering_agreement"))]
    with_jac = [base[0], ("Cell-type\nJaccard (mean)", g("celltype_jaccard_mean")), base[1], base[2]]

    for suffix, sc in [("", base), ("_jaccard", with_jac)]:
        fig, axes = plt.subplots(2, 2, figsize=(11, 9))

        ax = axes[0, 0]
        for i, n in enumerate(names):
            r = barcode_rank(M[n]); r = r[r > 0]
            ax.loglog(np.arange(1, len(r) + 1), r, label=f"{n} (n={len(r):,})", color=tool_color(n, i), lw=1.8)
        ax.set_xlabel("Barcode rank"); ax.set_ylabel("Total UMI per barcode")
        ax.set_title("Barcode rank"); ax.legend(loc="lower left")
        ax.grid(True, which="major", alpha=0.25, lw=0.5); panel(ax, "A")

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

        scorecard(axes[1, 0], [l for l, _ in sc], [v for _, v in sc], shared_n); panel(axes[1, 0], "C")

        ax = axes[1, 1]
        cts = list(next(iter(fr.values())).keys()) if fr else []
        x = np.arange(len(cts)); w = 0.8 / max(len(names), 1)
        for i, n in enumerate(names):
            ax.bar(x + i * w - 0.4 + w / 2, [fr.get(n, {}).get(c, 0) for c in cts], w, label=n, color=tool_color(n, i))
        ax.set_xticks(x); ax.set_xticklabels(cts, rotation=35, ha="right")
        ax.set_ylabel("Fraction of cells"); ax.set_title("Cell-type composition")
        ax.legend(loc="upper right"); panel(ax, "D")

        fig.tight_layout()
        print("saved", save(fig, os.path.join(outdir, "downstream_concordance" + suffix)), "(+ .pdf)")

if __name__ == "__main__":
    main()
