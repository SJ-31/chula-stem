#!/usr/bin/env ipython


import numpy as np
import pandas as pd

file = (
    "/home/shannc/Bio_SDD/chula-stem/analyses/output/HCC_abundance/hcc_expr_entrez.csv"
)


# * CUP-AI-Dx
# <2025-01-29 Wed> Finished the main function for prediction,
# just need to test out the classification performance
data_dir = "/home/shannc/Bio_SDD/chula-stem/extras/CUP-AI-Dx/data"

target_genes = pd.read_csv(
    f"{data_dir}/features_791.csv", header=None, index_col=0
).index

df = pd.read_csv(file, index_col=0)
if "X" not in str(df.index[0]):
    df.index = [f"X{i}" for i in df.index]

rows = []
for t in target_genes:
    try:
        rows.append(df.loc[[t],])
    except KeyError:
        rows.append(pd.DataFrame(0, columns=df.columns, index=[t]))
selected = pd.concat(rows, axis=0)


samples = list(df.columns)
transposed = np.transpose(df)

# * BPformer
sample_file = "/home/shannc/Bio_SDD/chula-stem/extras/BPformer/RNAseq/Raw/GEO-RNA-m.pkl"
sample = pd.read_pickle(sample_file)
