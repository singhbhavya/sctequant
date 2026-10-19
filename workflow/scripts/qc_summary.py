"""Collect STARsolo Summary.csv files into one table (one row per sample)."""

import os
import sys

import pandas as pd

out, *summaries = sys.argv[1:]
rows = {}
for path in summaries:
    sample = path.split(os.sep)[-5]  # <out>/<sample>/starsolo/Solo.out/Gene/Summary.csv
    rows[sample] = pd.read_csv(path, header=None, index_col=0).iloc[:, 0]
pd.DataFrame(rows).T.to_csv(out, index_label="sample")
