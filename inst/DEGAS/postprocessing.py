import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors


def centerFunc(x):
    """
    Center vector or matrix column-wise.
    """

    x = np.asarray(x, dtype=np.float32)

    if x.ndim == 1:
        return x - np.nanmean(x)

    return x - np.nanmean(x, axis=0, keepdims=True)


def toCorrCoeff(probs):
    """
    Convert DEGAS output probabilities [0, 1] to association scores [-1, 1].

    This follows the paper Eq.12.

    For binary classification:
        association = 2 * p - 1
    """

    arr = np.asarray(probs, dtype=np.float32)

    if arr.ndim == 1:
        k = 2
    else:
        k = arr.shape[1]
        if k < 2:
            k = 2

    return 2.0 * ((arr - 1.0 / k) / (2.0 - 2.0 / k) + 0.5) - 1.0


def knnSmooth(probs, locs, k=5):
    """
    kNN smoothing for DEGAS probability or association matrix.

    Parameters
    ----------
    probs : array-like
        n_cells × n_labels or n_cells vector.
    locs : array-like
        n_cells × n_dim coordinates, e.g. UMAP/tSNE.
    k : int
        Number of neighbors including itself.
    """

    probs_is_df = isinstance(probs, pd.DataFrame)
    probs_index = probs.index if probs_is_df else None
    probs_cols = probs.columns if probs_is_df else None

    p = np.asarray(probs, dtype=np.float32)
    locs = np.asarray(locs, dtype=np.float32)

    if p.ndim == 1:
        p2 = p.reshape(-1, 1)
    else:
        p2 = p

    k = int(k)
    k = max(1, min(k, locs.shape[0]))

    nn = NearestNeighbors(n_neighbors=k)
    nn.fit(locs)
    idx = nn.kneighbors(locs, return_distance=False)

    out = np.zeros_like(p2, dtype=np.float32)

    for i in range(p2.shape[0]):
        out[i, :] = np.nanmean(p2[idx[i], :], axis=0)

    if p.ndim == 1:
        out = out[:, 0]

    if probs_is_df:
        return pd.DataFrame(out, index=probs_index, columns=probs_cols)

    return out
