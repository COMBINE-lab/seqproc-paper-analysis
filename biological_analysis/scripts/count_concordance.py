#!/usr/bin/env python3
"""Count-level concordance across N tools: barcode-rank (knee), per-barcode UMI and per-gene
total Pearson correlation computed on log1p values (not raw counts). Reports every pairwise r.

  count_concordance.py <outdir> <name1>:<Gene_dir1> <name2>:<Gene_dir2> [<name3>:<dir3> ...]
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
import matplotlib; matplotlib.use("Agg")
import numpy as np, matplotlib.pyplot as plt
from numpy import log1p, corrcoef
from scipy.stats import spearmanr
from concordance_helpers import load_star_raw, per_barcode_umi, barcode_rank
from paper_style import set_paper_style, tool_color, panel, save

def main():
    set_paper_style()
    outdir = sys.argv[1]
    tools = [s.split(":", 1) for s in sys.argv[2:]]       # [(name, gene_dir), ...]
    os.makedirs(outdir, exist_ok=True)
    names = [n for n, _ in tools]
    M, BC, G = {}, {}, {}
    for n, d in tools:
        M[n], BC[n], G[n] = load_star_raw(d)
    ref = names[0]

    pbu = {n: dict(zip(BC[n], per_barcode_umi(M[n]))) for n in names}
    assert all(list(G[n]) == list(G[ref]) for n in names), "gene order differs across tools"
    pgt = {n: np.asarray(M[n].sum(0)).ravel() for n in names}

    def pair_r(a, b, per_barcode):
        if per_barcode:
            keys = sorted(set(pbu[a]) | set(pbu[b]))
            x = np.array([pbu[a].get(k, 0) for k in keys], float)
            y = np.array([pbu[b].get(k, 0) for k in keys], float)
        else:
            x, y = pgt[a], pgt[b]
        m = (x + y) > 0
        if m.sum() <= 1:
            return float("nan"), float("nan")
        pear = float(corrcoef(log1p(x[m]), log1p(y[m]))[0, 1])   # Pearson on log1p (dynamic range)
        spear = float(spearmanr(x[m], y[m]).correlation)         # Spearman (rank; transform-invariant)
        return pear, spear

    pairs = [(a, b) for i, a in enumerate(names) for b in names[i + 1:]]
    pb = {f"{a}|{b}": pair_r(a, b, True) for a, b in pairs}
    pg = {f"{a}|{b}": pair_r(a, b, False) for a, b in pairs}
    res = {
        "per_barcode_umi_pearson_logspace": {k: round(v[0], 4) for k, v in pb.items()},
        "per_barcode_umi_spearman":        {k: round(v[1], 4) for k, v in pb.items()},
        "per_gene_total_pearson_logspace": {k: round(v[0], 4) for k, v in pg.items()},
        "per_gene_total_spearman":         {k: round(v[1], 4) for k, v in pg.items()},
    }
    json.dump(res, open(os.path.join(outdir, "count_concordance.json"), "w"), indent=2)
    print("per-barcode  pearson(log):", res["per_barcode_umi_pearson_logspace"], "| spearman:", res["per_barcode_umi_spearman"])
    print("per-gene     pearson(log):", res["per_gene_total_pearson_logspace"], "| spearman:", res["per_gene_total_spearman"])

    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))

    # (A) barcode-rank knee, all tools
    for i, n in enumerate(names):
        r = barcode_rank(M[n]); r = r[r > 0]
        ax[0].loglog(np.arange(1, len(r) + 1), r, label=f"{n} (n={len(r):,})",
                     color=tool_color(n, i), lw=1.8)
    ax[0].set_xlabel("Barcode rank"); ax[0].set_ylabel("Total UMI per barcode")
    ax[0].set_title("Barcode rank"); ax[0].legend(loc="lower left")
    ax[0].grid(True, which="major", alpha=0.25, lw=0.5); panel(ax[0], "A")

    def scatter_panel(a, getx, gety, label, key, ttl):
        xr = getx(ref)
        for i, n in enumerate(names[1:], 1):
            y = gety(n); x = xr; m = (x + y) > 0
            pr = res[key][f"{ref}|{n}"]; sp = res[key.replace("pearson_logspace", "spearman")][f"{ref}|{n}"]
            a.scatter(log1p(x[m]), log1p(y[m]), s=6, alpha=0.35, color=tool_color(n, i),
                      edgecolors="none", label=f"{n}:  r={pr:.3f}, ρ={sp:.3f}")
        lim = max(log1p(xr).max(), 1); a.plot([0, lim], [0, lim], color="0.4", ls="--", lw=1)
        a.set_xlabel(f"seqproc  log1p({label})"); a.set_ylabel(f"other tool  log1p({label})")
        a.set_title(ttl); a.legend(loc="upper left", handletextpad=0.2)
        a.set_aspect("equal", adjustable="box")

    keys = sorted(set().union(*[set(pbu[n]) for n in names]))
    pbu_arr = {n: np.array([pbu[n].get(k, 0) for k in keys], float) for n in names}
    scatter_panel(ax[1], lambda n: pbu_arr[n], lambda n: pbu_arr[n],
                  "UMI / barcode", "per_barcode_umi_pearson_logspace", "Per-barcode UMI"); panel(ax[1], "B")
    scatter_panel(ax[2], lambda n: pgt[n], lambda n: pgt[n],
                  "gene UMI", "per_gene_total_pearson_logspace", "Per-gene total"); panel(ax[2], "C")

    fig.tight_layout()
    print("saved", save(fig, os.path.join(outdir, "count_concordance")), "(+ .pdf)")

if __name__ == "__main__":
    main()
