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
RANDOM_SEED = 0
np.random.seed(RANDOM_SEED)

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
        if present:
            sc.tl.score_genes(
                a, present, score_name=f"score_{ct}", random_state=RANDOM_SEED
            )
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
    sc.pp.pca(
        a,
        n_comps=min(n_pcs, a.n_obs - 1, a.n_vars - 1),
        random_state=RANDOM_SEED,
    )
    sc.pp.neighbors(a, n_neighbors=15, random_state=RANDOM_SEED)
    sc.tl.leiden(
        a,
        resolution=1.0,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        random_state=RANDOM_SEED,
    )
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
    sc.pp.pca(
        comb,
        n_comps=min(n_pcs, comb.n_obs - 1, comb.n_vars - 1),
        random_state=RANDOM_SEED,
    )
    sc.pp.neighbors(comb, n_neighbors=15, random_state=RANDOM_SEED)
    sc.tl.leiden(
        comb,
        resolution=1.0,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        random_state=RANDOM_SEED,
    )
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
    # pairwise cell-type agreement: fraction of shared cells each tool pair labels the same
    ct_agree_pw = {f"{a}|{b}": round(float(np.mean([type_lab[a][c] == type_lab[b][c] for c in shared])), 4) for a, b in pairs} if shared else {}

    # Per-cell-type Jaccard for every tool pair (label-set overlap; no
    # embedding).  Keep the pairwise values as first-class results: averaging
    # them is useful for a compact scorecard, but can conceal which comparison
    # drives a lower score.
    ct_jac = {}
    ct_jac_pw = {f"{a}|{b}": {} for a, b in pairs}
    for ct in MARKERS:
        js = []
        for a, b in pairs:
            sa = {c for c in shared if type_lab[a][c] == ct}
            sb = {c for c in shared if type_lab[b][c] == ct}
            u = sa | sb
            value = len(sa & sb) / len(u) if u else None
            ct_jac_pw[f"{a}|{b}"][ct] = round(float(value), 4) if value is not None else None
            if value is not None:
                js.append(value)
        ct_jac[ct] = round(float(np.mean(js)), 4) if js else None
    jvals = [v for v in ct_jac.values() if v is not None]
    mean_ct_jac = round(float(np.mean(jvals)), 4) if jvals else float("nan")
    ct_jac_pw_mean = {
        pair: round(float(np.mean([value for value in values.values() if value is not None])), 4)
        for pair, values in ct_jac_pw.items()
    }

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
        "tools": names, "min_umi": min_umi, "random_seed": RANDOM_SEED,
        "cells": {n: int(proc[n].n_obs) for n in names},
        "shared_cells": len(shared),
        "celltype_agreement_shared": round(type_agree, 4),
        "celltype_agreement_pairwise": ct_agree_pw,
        "celltype_jaccard_per_type": ct_jac,
        "celltype_jaccard_pairwise_per_type": ct_jac_pw,
        "celltype_jaccard_pairwise_mean": ct_jac_pw_mean,
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

    # ---- Jaccard supplement table (per cell type, for the supplement) ----
    with open(os.path.join(outdir, "jaccard_supplement.md"), "w") as fh:
        pair_keys = [f"{a}|{b}" for a, b in pairs]
        fh.write("| Cell type | " + " | ".join(pair.replace("|", " / ") for pair in pair_keys) + " | Mean |\n")
        fh.write("|---|" + "---:|" * (len(pair_keys) + 1) + "\n")
        for ct in MARKERS:
            values = [ct_jac_pw[pair].get(ct) for pair in pair_keys]
            rendered = [f"{value:.3f}" if value is not None else "n/a" for value in values]
            fh.write(f"| {ct} | " + " | ".join(rendered) + f" | {ct_jac[ct]:.3f} |\n")
        fh.write("| **Mean** | " + " | ".join(f"**{ct_jac_pw_mean[pair]:.3f}**" for pair in pair_keys) + f" | **{mean_ct_jac:.3f}** |\n")
    print("saved", os.path.join(outdir, "jaccard_supplement.md"))

    # ---- figure: cell-type composition + concordance scorecard, two variants ----
    def pairwise_bars(ax, include_jaccard):
        metrics = [
            ("Cell-type\nagreement", ct_agree_pw),
            ("Cluster ARI", ari),
        ]
        if include_jaccard:
            metrics.insert(1, ("Mean type\nJaccard", ct_jac_pw_mean))
        pair_keys = [f"{a}|{b}" for a, b in pairs]
        pair_labels = [pair.replace("|", " / ") for pair in pair_keys]
        pair_colors = ("#0072B2", "#D55E00", "#009E73")
        yb = np.arange(len(metrics))[::-1]
        height = 0.22
        for pair_index, (pair, pair_label, color) in enumerate(zip(pair_keys, pair_labels, pair_colors)):
            offset = (pair_index - (len(pair_keys) - 1) / 2) * height
            values = [metric_values[pair] for _, metric_values in metrics]
            bars = ax.barh(yb + offset, values, height, color=color, label=pair_label)
            for bar, value in zip(bars, values):
                ax.text(
                    min(value + 0.012, 1.005),
                    bar.get_y() + bar.get_height() / 2,
                    f"{value:.3f}",
                    va="center",
                    ha="left",
                    fontsize=8,
                )
        ax.axvline(1.0, color="0.5", ls="--", lw=1)
        ax.set_yticks(yb)
        ax.set_yticklabels([label for label, _ in metrics])
        # Agreement scores are bounded proportions; retain the zero baseline so
        # modest pairwise differences are not visually exaggerated.
        ax.set_xlim(0.0, 1.04)
        ax.set_xlabel("Pairwise score (1 = identical)")
        ax.set_title(f"Pairwise concordance\nshared cells n={len(shared)}", pad=36)
        ax.legend(
            loc="lower center",
            bbox_to_anchor=(0.5, 1.01),
            ncol=3,
            frameon=False,
            borderaxespad=0,
        )
        ax.spines["left"].set_visible(False)

    for suffix, include_jaccard in [("", False), ("_jaccard", True)]:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
        ax = axes[0]
        cts = list(MARKERS.keys()); x = np.arange(len(cts)); w = 0.8 / len(names)
        for i, n in enumerate(names):
            ax.bar(x + i * w - 0.4 + w / 2, [fracs[n][c] for c in cts], w, label=n, color=tool_color(n, i))
        ax.set_xticks(x); ax.set_xticklabels(cts, rotation=35, ha="right")
        ax.set_ylabel("Fraction of cells"); ax.set_title("Cell-type composition")
        ax.legend(loc="upper right"); panel(ax, "A")
        pairwise_bars(axes[1], include_jaccard); panel(axes[1], "B")
        fig.tight_layout()
        print("saved", save(fig, os.path.join(outdir, "biological_analysis" + suffix)), "(+ .pdf)")

if __name__ == "__main__":
    main()
