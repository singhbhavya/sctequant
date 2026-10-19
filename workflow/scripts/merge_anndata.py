"""Merge the per-sample AnnData files from one TE tool into one file."""

import sys

import anndata as ad
import pandas as pd

out, *inputs = sys.argv[1:]
adatas = [ad.read_h5ad(f) for f in inputs]
var = pd.concat([a.var for a in adatas])
var = var[~var.index.duplicated()]

merged = ad.concat(adatas, join="outer", fill_value=0, merge="first")
merged.var = var.loc[merged.var_names]
merged.uns["sctequant"] = adatas[0].uns["sctequant"]
for col in merged.obs.columns:
    if merged.obs[col].dtype == object:
        merged.obs[col] = merged.obs[col].astype("category")
merged.write_h5ad(out, compression="gzip")
print(f"Merged {len(inputs)} samples: {merged.n_obs} cells x {merged.n_vars} features")
