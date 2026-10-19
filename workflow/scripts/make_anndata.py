"""Combine STARsolo gene counts and one TE tool's counts into a single AnnData.

Cells are the cells STARsolo called (Gene/filtered). TE counts are matched to
those barcodes; cells with no TE reads get zeros.
"""

import argparse
import gzip
import os

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def open_maybe_gz(path):
    if os.path.exists(path):
        return path
    if os.path.exists(path + ".gz"):
        return path + ".gz"
    raise FileNotFoundError(path)


def read_lines(path):
    path = open_maybe_gz(path)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return [line.rstrip("\n").split("\t") for line in fh if line.strip()]


def read_mex(folder, matrix="matrix.mtx", features="features.tsv", barcodes="barcodes.tsv"):
    """Read a 10x-style matrix folder. Returns (cells x features csr, features rows, barcodes)."""
    mtx = scipy.io.mmread(open_maybe_gz(os.path.join(folder, matrix)))
    feats = read_lines(os.path.join(folder, features))
    bcs = [row[0] for row in read_lines(os.path.join(folder, barcodes))]
    return sp.csr_matrix(mtx.T, dtype=np.float32), feats, bcs


def load_solote(path, sample, level):
    folder = os.path.join(path, f"{sample}_{level}_MATRIX")
    X, feats, bcs = read_mex(folder)
    names = [row[0] for row in feats]
    keep = [i for i, n in enumerate(names) if n.startswith("SoloTE|")]
    return X[:, keep], [names[i][len("SoloTE|"):] for i in keep], bcs


def load_scte(path, te_names):
    scte = ad.read_h5ad(path)
    keep = np.array([n in te_names for n in scte.var_names])
    X = sp.csr_matrix(scte.X[:, keep], dtype=np.float32)
    return X, list(scte.var_names[keep]), list(scte.obs_names)


def load_irescue(path):
    X, feats, bcs = read_mex(path)
    return X, [row[0] for row in feats], bcs


def align_cells(X, barcodes, cells):
    """Reorder rows of X to `cells`, filling cells missing from X with zeros."""
    index = {bc: i for i, bc in enumerate(barcodes)}
    rows = [index.get(bc, -1) for bc in cells]
    present = [i for i, r in enumerate(rows) if r >= 0]
    out = sp.lil_matrix((len(cells), X.shape[1]), dtype=np.float32)
    if present:
        out[present, :] = X[[rows[i] for i in present], :]
    return out.tocsr()


def te_var(names, families, level):
    var = pd.DataFrame(index=pd.Index(names, name=None))
    var["gene_ids"] = names
    var["feature_type"] = "TE"
    if level == "locustes":
        # SoloTE locus names look like chr|start|end|subfamily:family:class|div|strand
        parts = [n.split("|") for n in names]
        sub = [p[3].split(":")[0] if len(p) > 3 else n for p, n in zip(parts, names)]
        var["te_subfamily"] = sub
    elif level in ("familytes", "classtes"):
        var["te_subfamily"] = np.nan
    else:
        var["te_subfamily"] = names
    var["te_family"] = var["te_subfamily"].map(families["family"])
    var["te_class"] = var["te_subfamily"].map(families["class"])
    if level == "familytes":
        var["te_family"] = names
    if level == "classtes":
        var["te_class"] = names
    return var


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", required=True)
    p.add_argument("--tool", required=True, choices=["solote", "scte", "irescue"])
    p.add_argument("--genes", required=True, help="STARsolo Solo.out/Gene/filtered folder")
    p.add_argument("--te", required=True, help="TE tool output")
    p.add_argument("--families", required=True)
    p.add_argument("--sample-sheet-row", required=True, help="CSV with one row of sample metadata")
    p.add_argument("--solote-level", default="subfamilytes")
    p.add_argument("--chemistry", default="")
    p.add_argument("--genome", default="")
    p.add_argument("--out", required=True)
    a = p.parse_args()

    families = pd.read_csv(a.families, sep="\t", index_col=0)

    G, gfeats, cells = read_mex(a.genes)
    gene_var = pd.DataFrame(
        {
            "gene_ids": [f[0] for f in gfeats],
            "feature_type": "Gene Expression",
            "te_subfamily": np.nan,
            "te_family": np.nan,
            "te_class": np.nan,
        },
        index=pd.Index([f[1] if len(f) > 1 else f[0] for f in gfeats]),
    )

    level = a.solote_level
    if a.tool == "solote":
        T, tnames, tbcs = load_solote(a.te, a.sample, level)
    elif a.tool == "scte":
        T, tnames, tbcs = load_scte(a.te, set(families.index))
        level = "subfamilytes"
    else:
        T, tnames, tbcs = load_irescue(a.te)
        level = "subfamilytes"

    T = align_cells(T, tbcs, cells)
    var = pd.concat([gene_var, te_var(tnames, families, level)])
    X = sp.hstack([G, T], format="csr", dtype=np.float32)

    meta = pd.read_csv(a.sample_sheet_row, dtype=str).iloc[0].to_dict()
    obs = pd.DataFrame(index=pd.Index([f"{a.sample}_{bc}" for bc in cells]))
    obs["barcode"] = cells
    for key, value in meta.items():
        obs[key] = value
    obs["sample"] = a.sample
    obs["gene_counts"] = np.asarray(G.sum(axis=1)).ravel()
    obs["te_counts"] = np.asarray(T.sum(axis=1)).ravel()
    total = obs["gene_counts"] + obs["te_counts"]
    obs["pct_te"] = np.where(total > 0, 100 * obs["te_counts"] / total.where(total > 0, 1), 0.0)
    for col in obs.columns:
        if obs[col].dtype == object:
            obs[col] = obs[col].astype("category")

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names_make_unique()
    for col in ("feature_type", "te_family", "te_class"):
        adata.var[col] = adata.var[col].astype("category")
    adata.uns["sctequant"] = {
        "te_tool": a.tool,
        "te_level": level.replace("tes", ""),
        "chemistry": a.chemistry,
        "genome": a.genome,
    }
    adata.write_h5ad(a.out, compression="gzip")

    print(
        f"{a.sample} [{a.tool}]: {adata.n_obs} cells, {G.shape[1]} genes, {T.shape[1]} TE features, "
        f"median {np.median(obs['pct_te']):.2f}% TE UMIs per cell"
    )


if __name__ == "__main__":
    main()
