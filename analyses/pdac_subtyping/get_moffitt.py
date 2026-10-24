#!/usr/bin/env ipython


import anndata as ad
import scanpy as sc
from chula_stem import pdac_subtyping as pds
from pyhere import here

adatas = []
for file in ["crc-pdac.h5ad", "crc-pdac_test.h5ad"]:
    tmp = ad.read_h5ad(
        here("analyses", "data_all", "repos", "too-predict", "training", file)
    )
    tmp = tmp[
        (tmp.obs["Project_ID"].str.contains("CHULA"))
        & (tmp.obs["tumor_type"] == "PDAC"),
        :,
    ]
    adatas.append(tmp)

adata = ad.concat(adatas, axis="obs", join="inner", merge="first")
workdir = here("analyses", "pdac_subtyping")
adata.write_h5ad(workdir / "cohort.h5ad")

model = pds.MoffittBasal(threshold=0.44)  # From running last time
model.fit(adata=adata, gene_col="GENENAME")
model.predict(inplace=True)

adata.obs.loc[:, ["Case_ID", "moffitt_is_basal"]].rename(
    {"moffitt_is_basal": "label"}, axis="columns"
).to_csv(workdir / "moffitt_is_basal.tsv", index=False, sep="\t", na_rep="NA")
