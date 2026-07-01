#!/usr/bin/env python3
"""Phase 2A biological readout, N-tool. Cell-calling, clustering, cell typing, and tool
concordance reported ENTIRELY as metrics (no UMAP/t-SNE visualization, which requires visual
judgment and is contested). Concordance is quantified by cell-type agreement, per-cell-type
Jaccard, between-tool clustering ARI, and joint co-clustering, all computed from counts and
kNN-graph Leiden clusterings (never from a 2D embedding).

  biological_analysis.py <outdir> <min_umi> <name1>:<Gene_dir1> <name2>:<Gene_dir2> [<name3>:<dir3> ...]
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

MARKERS = {
    "Neuron":      ["Snap25", "Syt1", "Rbfox3", "Meg3"],
    "Astrocyte":   ["Gfap", "Aqp4", "Slc1a3", "Aldoc"],
    "Oligodendro": ["Plp1", "Mbp", "Mog", "Mobp"],
    "OPC":         ["Pdgfra", "Cspg4"],
    "Microglia":   ["C1qa", "C1qb", "Cx3cr1", "Csf1r"],
    "Endothelial": ["Cldn5", "Flt1", "Pecam1"],
}

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
    type_cells(a)
    return a

def joint_coclustering(raw, shared, min_umi, n_hvg=2000, n_pcs=30):
    """Embed shared cells from all tools in ONE kNN graph, cluster with Leiden (no UMAP), and
    return (co_clustering_agreement, tool_mixing_entropy). Co-clustering = fraction of cells
    whose per-tool versions land in the same joint cluster; mixing = mean normalized per-cluster
    tool entropy (1 = perfectly mixed)."""
    mats = []
    for tool, a in raw.items():
        sub = a[shared].copy()
        sub.obs["tool"] = tool
        sub.obs["cell"] = list(shared)
        sub.obs_names = [f"{tool}|{b}" for b in shared]
        mats.append(sub)
    comb = ad.concat(mats)
    sc.pp.normalize_total(comb, target_sum=1e4); sc.pp.log1p(comb)
    sc.pp.highly_variable_genes(comb, n_top_genes=min(n_hvg, comb.n_vars - 1))
    sc.pp.pca(comb, n_comps=min(n_pcs, comb.n_obs - 1, comb.n_vars - 1))
    sc.pp.neighbors(comb, n_neighbors=15)
    sc.tl.leiden(comb, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
    import collections
    by_cell = collections.defaultdict(list)
    for cl, cell in zip(comb.obs["leiden"], comb.obs["cell"]):
        by_cell[cell].append(cl)
    co = float(np.mean([len(set(v)) == 1 for v in by_cell.values()]))
    tools = sorted(set(comb.obs["tool"]))
    ent = []
    for cl in comb.obs["leiden"].unique():
        sub = comb.obs[comb.obs["leiden"] == cl]["tool"]
        p = np.array([np.mean(sub == t) for t in tools]); p = p[p > 0]
        h = -(p * np.log(p)).sum() / np.log(len(tools)) if len(tools) > 1 else 1.0
        ent.append((len(sub), h))
    w = np.array([n for n, _ in ent], float)
    mixing = float(np.average([h for _, h in ent], weights=w))
    return co, mixing

def main():
    outdir, min_umi = sys.argv[1], int(sys.argv[2])
    tools = [s.split(":", 1) for s in sys.argv[3:]]                 # [(name, gene_dir), ...]
    os.makedirs(outdir, exist_ok=True)
    raw = {n: load(d) for n, d in tools}                            # for joint co-clustering
    proc = {n: process(load(d), min_umi=min_umi) for n, d in tools} # per-tool clustering/typing
    names = [n for n, _ in tools]
    for n in names:
        print(f"{n}: {proc[n].n_obs} cells, {proc[n].n_vars} genes, {proc[n].obs['leiden'].nunique()} clusters")

    called = {n: set(proc[n].obs_names) for n in names}
    shared = sorted(set.intersection(*called.values()))
    fracs = {n: {ct: float((proc[n].obs["cell_type"] == ct).mean()) for ct in MARKERS} for n in names}
    type_lab = {n: dict(zip(proc[n].obs_names, proc[n].obs["cell_type"])) for n in names}
    ref = names[0]
    pairs = [(a, b) for i, a in enumerate(names) for b in names[i + 1:]]

    # cell-type agreement: fraction of shared cells with the same label in EVERY tool
    type_agree = float(np.mean([all(type_lab[n][c] == type_lab[ref][c] for n in names) for c in shared])) if shared else float("nan")

    # per-cell-type Jaccard across tool pairs (label-set overlap; no embedding)
    ct_jac = {}
    for ct in MARKERS:
        js = []
        for a, b in pairs:
            sa = {c for c in shared if type_lab[a][c] == ct}
            sb = {c for c in shared if type_lab[b][c] == ct}
            u = sa | sb
            if u: js.append(len(sa & sb) / len(u))
        ct_jac[ct] = round(float(np.mean(js)), 4) if js else None
    jvals = [v for v in ct_jac.values() if v is not None]
    mean_ct_jac = round(float(np.mean(jvals)), 4) if jvals else float("nan")

    # between-tool clustering ARI (independent Leiden; full-depth meaningful, low-depth fragile)
    ari = {}
    if len(shared) > 10:
        leid = {n: proc[n][shared].obs["leiden"].astype(str).values for n in names}
        for a, b in pairs:
            ari[f"{a}|{b}"] = round(float(adjusted_rand_score(leid[a], leid[b])), 4)
    mean_ari = round(float(np.mean(list(ari.values()))), 4) if ari else float("nan")

    # joint co-clustering (kNN-graph Leiden, no UMAP)
    co_cluster, mixing = (float("nan"), float("nan"))
    if len(shared) > 10:
        co_cluster, mixing = joint_coclustering(raw, shared, min_umi)

    metrics = {
        "tools": names, "min_umi": min_umi,
        "cells": {n: int(proc[n].n_obs) for n in names},
        "shared_cells": len(shared),
        "celltype_agreement_shared": round(type_agree, 4),
        "celltype_jaccard_per_type": ct_jac,
        "celltype_jaccard_mean": mean_ct_jac,
        "celltype_fractions": {n: {k: round(v, 4) for k, v in fracs[n].items()} for n in names},
        "cluster_ari_pairwise": ari,
        "cluster_ari_mean": mean_ari,
        "joint_coclustering_agreement": round(co_cluster, 4),
        "tool_mixing_entropy": round(mixing, 4),
    }
    json.dump(metrics, open(os.path.join(outdir, "biological_metrics.json"), "w"), indent=2)
    print("celltype_agree", metrics["celltype_agreement_shared"],
          "| celltype_jaccard(mean)", mean_ct_jac,
          "| cluster ARI(mean)", mean_ari,
          "| co-clustering", metrics["joint_coclustering_agreement"],
          "| tool-mixing", metrics["tool_mixing_entropy"])

    # ---- figure: cell-type composition + quantitative concordance scorecard (no embedding) ----
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))

    ax = axes[0]
    cts = list(MARKERS.keys()); x = np.arange(len(cts)); w = 0.8 / len(names)
    for i, n in enumerate(names):
        ax.bar(x + i * w - 0.4 + w / 2, [fracs[n][c] for c in cts], w, label=n, color=tool_color(n, i))
    ax.set_xticks(x); ax.set_xticklabels(cts, rotation=35, ha="right")
    ax.set_ylabel("Fraction of cells"); ax.set_title("Cell-type composition")
    ax.legend(loc="upper right"); panel(ax, "A")

    ax = axes[1]
    labels = ["Cell-type\nagreement", "Cell-type\nJaccard (mean)", "Cluster ARI\n(mean)", "Co-clustering\n(joint)"]
    vals = [type_agree, mean_ct_jac, mean_ari, co_cluster]
    yb = np.arange(len(labels))[::-1]
    ax.barh(yb, vals, color="#1b7837", height=0.55)
    ax.axvline(1.0, color="0.5", ls="--", lw=1)
    for y_, v in zip(yb, vals):
        if v == v:
            ax.text(min(v, 0.97) - 0.02, y_, f"{v:.3f}", va="center", ha="right", color="white", fontweight="bold")
    ax.set_yticks(yb); ax.set_yticklabels(labels); ax.set_xlim(0, 1.08); ax.set_xlabel("Score (1 = identical)")
    ax.set_title(f"Tool concordance\nshared cells n={len(shared)}, {len(names)} tools")
    ax.spines["left"].set_visible(False); panel(ax, "B")

    fig.tight_layout()
    print("saved", save(fig, os.path.join(outdir, "biological_analysis")), "(+ .pdf)")

if __name__ == "__main__":
    main()
