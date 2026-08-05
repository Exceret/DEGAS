import numpy as np


def normFunc(x):
    """
    DEGAS-style z-score for one vector.
    Kept for compatibility with original DEGAS naming.
    """

    x = np.asarray(x, dtype=np.float32)
    mu = np.nanmean(x)
    sd = np.nanstd(x)

    if not np.isfinite(sd) or sd == 0:
        sd = 1.0

    return (x - mu) / sd


def scaleFunc(x):
    """
    Scale one vector to [0, 1].
    Kept for compatibility with original DEGAS naming.
    """

    x = np.asarray(x, dtype=np.float32)
    xmin = np.nanmin(x)
    xmax = np.nanmax(x)
    denom = xmax - xmin

    if not np.isfinite(denom) or denom == 0:
        denom = 1.0

    return (x - xmin) / denom


def scale_expression(x):
    """
    Matrix version used internally.

    Input
    -----
    x : samples/cells × genes

    Steps
    -----
    1. sample-wise z-score
    2. sample-wise [0, 1] scaling
    """

    x = np.asarray(x, dtype=np.float32)
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)

    mean = x.mean(axis=1, keepdims=True)
    std = x.std(axis=1, keepdims=True)
    std[std == 0] = 1.0

    z = (x - mean) / std

    zmin = z.min(axis=1, keepdims=True)
    zmax = z.max(axis=1, keepdims=True)
    denom = zmax - zmin
    denom[denom == 0] = 1.0

    return ((z - zmin) / denom).astype(np.float32)


def preprocessCounts(X):
    """
    Compatibility function with original DEGAS.

    Original DEGAS expects:
        rows = genes
        columns = cells/samples

    Returns:
        rows = cells/samples
        columns = genes

    This function applies:
        log2(X + 1)
        transpose
        sample-wise z-score
        sample-wise [0, 1] scaling
    """

    X = np.asarray(X, dtype=np.float32)
    X = np.log2(X + 1.0)
    X = X.T
    return scale_expression(X)


def align_genes(sc_matrix, bulk_matrix, sc_gene_names, bulk_gene_names):
    """
    Align scRNA and bulk matrices by common gene names.

    Input matrices:
        sc_matrix   : cells × genes
        bulk_matrix : samples × genes

    Returns
    -------
    sc_aligned
    bulk_aligned
    common_genes
    """

    sc_matrix = np.asarray(sc_matrix, dtype=np.float32)
    bulk_matrix = np.asarray(bulk_matrix, dtype=np.float32)

    sc_gene_names = np.asarray(sc_gene_names).astype(str)
    bulk_gene_names = np.asarray(bulk_gene_names).astype(str)

    if len(sc_gene_names) != sc_matrix.shape[1]:
        raise ValueError("len(sc_gene_names) must equal sc_matrix.shape[1].")

    if len(bulk_gene_names) != bulk_matrix.shape[1]:
        raise ValueError("len(bulk_gene_names) must equal bulk_matrix.shape[1].")

    sc_map = {}
    for i, g in enumerate(sc_gene_names):
        if g not in sc_map:
            sc_map[g] = i

    bulk_map = {}
    for i, g in enumerate(bulk_gene_names):
        if g not in bulk_map:
            bulk_map[g] = i

    common_genes = sorted(set(sc_map.keys()) & set(bulk_map.keys()))

    if len(common_genes) == 0:
        raise ValueError("No common genes found between scRNA and bulk matrices.")

    sc_idx = [sc_map[g] for g in common_genes]
    bulk_idx = [bulk_map[g] for g in common_genes]

    return sc_matrix[:, sc_idx], bulk_matrix[:, bulk_idx], common_genes


def rank_normalize(x):
    """
    Convert vector to [0, 1] rank score.
    Larger original value receives larger rank score.
    """

    x = np.asarray(x, dtype=np.float32)
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)

    if len(x) <= 1:
        return np.zeros_like(x, dtype=np.float32)

    order = np.argsort(x)
    ranks = np.empty_like(order, dtype=np.float32)
    ranks[order] = np.arange(len(x), dtype=np.float32)

    return ranks / float(len(x) - 1)
