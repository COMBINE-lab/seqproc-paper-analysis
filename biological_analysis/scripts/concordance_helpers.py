"""Load STARsolo raw matrices and compare two tools' outputs.
Dependencies: scipy.io, numpy only (no scanpy/pandas)."""
import gzip, os
import numpy as np
from scipy.io import mmread

def _open(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

def _pick(d, name):
    p = os.path.join(d, name)
    return p if os.path.exists(p) else p + ".gz"

def load_star_raw(solo_gene_dir):
    """Return (cells x genes CSR, barcodes ndarray, gene_names ndarray) from a STARsolo Gene dir."""
    raw = os.path.join(solo_gene_dir, "raw")
    M = mmread(_pick(raw, "matrix.mtx")).tocsr()        # genes x cells
    bc = np.array([l.strip() for l in _open(_pick(raw, "barcodes.tsv"))])
    feats = [l.rstrip("\n").split("\t") for l in _open(_pick(raw, "features.tsv"))]
    genes = np.array([f[1] if len(f) > 1 else f[0] for f in feats])
    return M.T.tocsr(), bc, genes                        # cells x genes

def per_barcode_umi(cxg):
    return np.asarray(cxg.sum(axis=1)).ravel()

def barcode_rank(cxg):
    return np.sort(per_barcode_umi(cxg))[::-1]

def aligned_barcode_umi(bc_a, umi_a, bc_b, umi_b):
    """Per-barcode UMI for tools a and b aligned on the union of barcodes."""
    da = dict(zip(bc_a, umi_a)); db = dict(zip(bc_b, umi_b))
    keys = sorted(set(da) | set(db))
    a = np.array([da.get(k, 0) for k in keys], float)
    b = np.array([db.get(k, 0) for k in keys], float)
    return np.array(keys), a, b
