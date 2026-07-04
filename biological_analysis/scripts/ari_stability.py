#!/usr/bin/env python3
"""Stability of the between-tool cluster ARI.

The main analysis reports one pairwise ARI per tool pair from a single Leiden run. That number is
sensitive to Leiden's random seed and to the clustering resolution, so on its own it cannot tell us
whether the splitcode-vs-others ARI gap is a reproducible property of the data or run-to-run noise.
This script answers that: it preprocesses each tool once (identically to biological_analysis.py),
then reclusters across many Leiden seeds (and, optionally, resolutions) and reports the spread of
the pairwise ARI. If the seqproc/matchbox pair and the splitcode pairs overlap across seeds, the gap
is within the metric's variability; if they stay separated, the gap is real and tracks recovery.

  ari_stability.py <outdir> <n_seeds> <min_umi> <name1>:<Gene_dir1> <name2>:<dir2> ... [--res=0.5,1.0,1.5]

Clustering matches the paper exactly: flavor="igraph", n_iterations=2, directed=False, resolution 1.0
by default. Writes ari_stability.json.
"""
import sys, os, json, gzip
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import scanpy as sc
import anndata as ad
from scipy.io import mmread
from sklearn.metrics import adjusted_rand_score

sc.settings.verbosity = 0

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

def preprocess(a, min_umi, n_hvg=2000, n_pcs=30):
    """Deterministic pipeline through the kNN graph; only Leiden below varies with the seed."""
    sc.pp.filter_cells(a, min_counts=min_umi)
    sc.pp.filter_genes(a, min_cells=3)
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=min(n_hvg, a.n_vars - 1))
    sc.pp.pca(a, n_comps=min(n_pcs, a.n_obs - 1, a.n_vars - 1))
    sc.pp.neighbors(a, n_neighbors=15)
    return a

def stats(v):
    v = np.asarray(v, float)
    return {"min": round(float(v.min()), 4), "median": round(float(np.median(v)), 4),
            "max": round(float(v.max()), 4), "mean": round(float(v.mean()), 4),
            "std": round(float(v.std()), 4), "n": int(v.size)}

def main():
    pos = [x for x in sys.argv[1:] if not x.startswith("--")]
    outdir, n_seeds, min_umi = pos[0], int(pos[1]), int(pos[2])
    tools = [s.split(":", 1) for s in pos[3:]]
    res_flag = [x for x in sys.argv[1:] if x.startswith("--res")]
    resolutions = [float(r) for r in res_flag[0].split("=", 1)[1].split(",")] if res_flag else [1.0]
    os.makedirs(outdir, exist_ok=True)
    names = [n for n, _ in tools]

    proc = {n: preprocess(load(d), min_umi) for n, d in tools}
    called = {n: set(proc[n].obs_names) for n in names}
    shared = sorted(set.intersection(*called.values()))
    pairs = [(a, b) for i, a in enumerate(names) for b in names[i + 1:]]
    print(f"shared cells: {len(shared)}; seeds: {n_seeds}; resolutions: {resolutions}")

    # seed-only spread at each resolution, keyed "<pair>@res<r>", plus the reference (seed 0, res 1.0)
    per = {f"{a}|{b}@res{r}": [] for a, b in pairs for r in resolutions}
    ref = {}
    for r in resolutions:
        for s in range(n_seeds):
            lab = {}
            for n in names:
                sc.tl.leiden(proc[n], resolution=r, random_state=s, flavor="igraph",
                             n_iterations=2, directed=False, key_added="l")
                lab[n] = proc[n][shared].obs["l"].astype(str).values
            for a, b in pairs:
                ari = float(adjusted_rand_score(lab[a], lab[b]))
                per[f"{a}|{b}@res{r}"].append(ari)
                if s == 0 and r == 1.0:
                    ref[f"{a}|{b}"] = round(ari, 4)

    summary = {
        "shared_cells": len(shared), "n_seeds": n_seeds, "resolutions": resolutions,
        "clustering": "leiden flavor=igraph n_iterations=2 directed=False",
        "reference_seed0_res1.0": ref,                               # should match the paper's table
        "seed_spread_per_resolution": {k: stats(v) for k, v in per.items()},
    }
    # combined spread over all seeds x resolutions, per pair
    summary["seed_x_resolution_spread"] = {
        f"{a}|{b}": stats([v for r in resolutions for v in per[f"{a}|{b}@res{r}"]]) for a, b in pairs
    }
    json.dump(summary, open(os.path.join(outdir, "ari_stability.json"), "w"), indent=2)

    print("reference (seed 0, res 1.0):", ref)
    print("seed spread at res 1.0:")
    for a, b in pairs:
        st = summary["seed_spread_per_resolution"][f"{a}|{b}@res1.0"]
        print(f"  {a}|{b}: median={st['median']}  [min {st['min']}, max {st['max']}]  mean={st['mean']}±{st['std']}")
    print("saved", os.path.join(outdir, "ari_stability.json"))

if __name__ == "__main__":
    main()
