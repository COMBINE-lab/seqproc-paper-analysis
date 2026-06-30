#!/usr/bin/env python3
"""Phase 2A biological readout, N-tool. Cell-calling, clustering, cell typing, and
tool concordance. Leads with the depth-robust metrics (cell-type fractions/agreement)
and uses a JOINT-EMBEDDING co-clustering metric instead of fragile independent-Leiden ARI.

  biological_analysis.py <outdir> <min_umi> <name1>:<Gene_dir1> <name2>:<Gene_dir2> [<name3>:<dir3> ...]

Joint-embedding concordance: the shared called cells are embedded ONCE in a combined
space (each cell contributed by every tool, tagged by tool), clustered jointly, then we
ask (a) do a cell's per-tool versions land in the same joint cluster (co-clustering
agreement), and (b) are clusters tool-mixed rather than tool-segregated (mixing entropy).
Independent-Leiden ARI is still reported but flagged as full-data-only / fragile.
"""
import sys, os, json, gzip
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
from scipy.io import mmread
from sklearn.metrics import adjusted_rand_score
from paper_style import set_paper_style, tool_color, panel, save

sc.settings.verbosity = 0
set_paper_style()

# fixed colors for the coarse cell types (consistent across panels)
CT_COLORS = {"Neuron": "#4C72B0", "Astrocyte": "#DD8452", "Oligodendro": "#55A868",
             "OPC": "#C44E52", "Microglia": "#8172B3", "Endothelial": "#937860", "Unknown": "#BBBBBB"}

MARKERS = {
    "Neuron":      ["Snap25", "Syt1", "Rbfox3", "Meg3"],
    "Astrocyte":   ["Gfap", "Aqp4", "Slc1a3", "Aldoc"],
    "Oligodendro": ["Plp1", "Mbp", "Mog", "Mobp"],
    "OPC":         ["Pdgfra", "Cspg4"],
    "Microglia":   ["C1qa", "C1qb", "Cx3cr1", "Csf1r"],
    "Endothelial": ["Cldn5", "Flt1", "Pecam1"],
}
PALETTE = ["#1565c0", "#e8743b", "#2e7d32", "#8e44ad", "#c0392b"]

def _pick(d, name):
    p = os.path.join(d, name); return p if os.path.exists(p) else p + ".gz"
def _open(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

def load(gene_dir):
    raw = os.path.join(gene_dir, "raw")
    M = mmread(_pick(raw, "matrix.mtx")).tocsr().T.tocsr()      # cells x genes
    bc = np.array([l.strip() for l in _open(_pick(raw, "barcodes.tsv"))])
    feats = [l.rstrip("\n").split("\t") for l in _open(_pick(raw, "features.tsv"))]
    genes = np.array([f[1] if len(f) > 1 else f[0] for f in feats])
    a = ad.AnnData(M.astype("float32")); a.obs_names = bc; a.var_names = genes
    a.var_names_make_unique()
    return a

def type_cells(a):
    for ct, mk in MARKERS.items():
        present = [g for g in mk if g in a.var_names]
        if present: sc.tl.score_genes(a, present, score_name=f"score_{ct}")
    cts = [ct for ct in MARKERS if f"score_{ct}" in a.obs]
    if cts:
        S = np.vstack([a.obs[f"score_{ct}"].values for ct in cts])
        a.obs["cell_type"] = np.array(cts)[S.argmax(0)]
    else:
        a.obs["cell_type"] = "Unknown"

def process(a, min_umi=200, n_hvg=2000, n_pcs=30):
    sc.pp.filter_cells(a, min_counts=min_umi)
    sc.pp.filter_genes(a, min_cells=3)
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=min(n_hvg, a.n_vars - 1))
    sc.pp.pca(a, n_comps=min(n_pcs, a.n_obs - 1, a.n_vars - 1))
    sc.pp.neighbors(a, n_neighbors=15)
    sc.tl.leiden(a, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(a)
    type_cells(a)
    return a

def joint_embedding(raw, shared, min_umi, n_hvg=2000, n_pcs=30):
    """Embed shared cells from all tools in one space; return combined AnnData + metrics."""
    mats = []
    for tool, a in raw.items():
        sub = a[shared].copy()
        sub.obs["tool"] = tool
        sub.obs["cell"] = list(shared)
        sub.obs_names = [f"{tool}|{b}" for b in shared]
        mats.append(sub)
    comb = ad.concat(mats)
    comb.layers["counts"] = comb.X.copy()
    sc.pp.normalize_total(comb, target_sum=1e4); sc.pp.log1p(comb)
    sc.pp.highly_variable_genes(comb, n_top_genes=min(n_hvg, comb.n_vars - 1))
    sc.pp.pca(comb, n_comps=min(n_pcs, comb.n_obs - 1, comb.n_vars - 1))
    sc.pp.neighbors(comb, n_neighbors=15)
    sc.tl.leiden(comb, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(comb)
    # (a) co-clustering: per cell, do all tools' versions share a joint cluster?
    import collections
    by_cell = collections.defaultdict(list)
    for cl, cell in zip(comb.obs["leiden"], comb.obs["cell"]):
        by_cell[cell].append(cl)
    co = np.mean([len(set(v)) == 1 for v in by_cell.values()])
    # (b) tool mixing: mean normalized tool-entropy per cluster (1 = perfectly mixed)
    tools = sorted(set(comb.obs["tool"]))
    ent = []
    for cl in comb.obs["leiden"].unique():
        sub = comb.obs[comb.obs["leiden"] == cl]["tool"]
        p = np.array([np.mean(sub == t) for t in tools]); p = p[p > 0]
        h = -(p * np.log(p)).sum() / np.log(len(tools)) if len(tools) > 1 else 1.0
        ent.append((len(sub), h))
    w = np.array([n for n, _ in ent], float)
    mixing = float(np.average([h for _, h in ent], weights=w))
    return comb, float(co), mixing

def main():
    outdir, min_umi = sys.argv[1], int(sys.argv[2])
    tools = [s.split(":", 1) for s in sys.argv[3:]]                 # [(name, gene_dir), ...]
    os.makedirs(outdir, exist_ok=True)
    raw = {n: load(d) for n, d in tools}                            # for joint embedding
    proc = {n: process(load(d), min_umi=min_umi) for n, d in tools} # per-tool clustering/typing
    names = [n for n, _ in tools]
    for n in names:
        print(f"{n}: {proc[n].n_obs} cells, {proc[n].n_vars} genes, {proc[n].obs['leiden'].nunique()} clusters")

    # shared called cells = passed min_umi in every tool
    called = {n: set(proc[n].obs_names) for n in names}
    shared = sorted(set.intersection(*called.values()))

    # depth-robust: cell-type fractions + agreement
    fracs = {n: {ct: float((proc[n].obs["cell_type"] == ct).mean()) for ct in MARKERS} for n in names}
    type_lab = {n: dict(zip(proc[n].obs_names, proc[n].obs["cell_type"])) for n in names}
    ref = names[0]
    type_agree = float(np.mean([all(type_lab[n][c] == type_lab[ref][c] for n in names) for c in shared])) if shared else float("nan")

    # joint-embedding concordance (robust clustering metric)
    co_cluster, mixing = (float("nan"), float("nan")); comb = None
    if len(shared) > 10:
        comb, co_cluster, mixing = joint_embedding(raw, shared, min_umi)

    # fragile independent-Leiden ARI, per tool-pair (kept, flagged; depth-sensitive)
    ari = {}
    if len(shared) > 10:
        leid = {n: proc[n][shared].obs["leiden"].astype(str).values for n in names}
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                ari[f"{a}|{b}"] = round(float(adjusted_rand_score(leid[a], leid[b])), 4)

    metrics = {
        "tools": names, "min_umi": min_umi,
        "cells": {n: int(proc[n].n_obs) for n in names},
        "shared_cells": len(shared),
        "depth_robust": {
            "celltype_agreement_shared": round(type_agree, 4),
            "celltype_fractions": {n: {k: round(v, 4) for k, v in fracs[n].items()} for n in names},
        },
        "clustering_joint_embedding": {
            "co_clustering_agreement": round(co_cluster, 4),   # primary, robust
            "tool_mixing_entropy": round(mixing, 4),           # 1.0 = perfectly mixed
        },
        "clustering_independent_ari_FRAGILE": ari,             # per-pair; depth-sensitive, full-data only
    }
    json.dump(metrics, open(os.path.join(outdir, "biological_metrics.json"), "w"), indent=2)
    print("metrics:", json.dumps(metrics["clustering_joint_embedding"]),
          "| celltype_agree", metrics["depth_robust"]["celltype_agreement_shared"],
          "| ARI(fragile)", metrics["clustering_independent_ari_FRAGILE"])

    # ---- paper figure: 2x2, depth-robust + joint embedding ----
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    # (A) cell-type fractions across tools
    ax = axes[0, 0]
    cts = list(MARKERS.keys()); x = np.arange(len(cts)); w = 0.8 / len(names)
    for i, n in enumerate(names):
        ax.bar(x + i * w - 0.4 + w / 2, [fracs[n][c] for c in cts], w, label=n, color=tool_color(n, i))
    ax.set_xticks(x); ax.set_xticklabels(cts, rotation=35, ha="right")
    ax.set_ylabel("Fraction of cells"); ax.set_title("Cell-type composition")
    ax.legend(loc="upper right"); panel(ax, "A")

    # (B) joint embedding colored by tool
    ax = axes[0, 1]
    if comb is not None:
        u = comb.obsm["X_umap"]
        for i, n in enumerate(names):
            m = comb.obs["tool"].values == n
            ax.scatter(u[m, 0], u[m, 1], s=10, alpha=0.7, label=n, color=tool_color(n, i), edgecolors="none")
        ax.legend(loc="best", markerscale=1.6)
    ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("Joint embedding by tool  (overlap = agreement)"); panel(ax, "B")

    # (C) joint embedding colored by cell type
    ax = axes[1, 0]
    if comb is not None:
        # transfer per-tool cell_type onto the joint cells (by tool|barcode)
        lab = {f"{n}|{b}": t for n in names for b, t in zip(proc[n].obs_names, proc[n].obs["cell_type"])}
        ctj = [lab.get(name, "Unknown") for name in comb.obs_names]
        for ct in MARKERS:
            m = np.array([c == ct for c in ctj])
            if m.any():
                ax.scatter(u[m, 0], u[m, 1], s=10, alpha=0.7, label=ct, color=CT_COLORS[ct], edgecolors="none")
        ax.legend(loc="best", markerscale=1.6, ncol=1)
    ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("Joint embedding by cell type"); panel(ax, "C")

    # (D) robust concordance metrics
    ax = axes[1, 1]
    labels = ["Cell-type\nagreement", "Co-clustering\n(joint)", "Tool-mixing\nentropy"]
    vals = [type_agree, co_cluster, mixing]
    yb = np.arange(len(labels))[::-1]
    ax.barh(yb, vals, color="#1b7837", height=0.55)
    ax.axvline(1.0, color="0.5", ls="--", lw=1)
    for y_, v in zip(yb, vals):
        ax.text(min(v, 0.98) - 0.02, y_, f"{v:.3f}", va="center", ha="right", color="white", fontweight="bold")
    ax.set_yticks(yb); ax.set_yticklabels(labels); ax.set_xlim(0, 1.08); ax.set_xlabel("Score (1 = identical)")
    ax.set_title(f"Tool concordance (robust)\nshared cells n={len(shared)}, {len(names)} tools")
    ax.spines["left"].set_visible(False); panel(ax, "D")

    fig.tight_layout()
    print("saved", save(fig, os.path.join(outdir, "biological_analysis")), "(+ .pdf)")

if __name__ == "__main__":
    main()
