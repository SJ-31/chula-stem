#!/usr/bin/env ipython
import re
from pathlib import Path

import anndata as ad
import chula_stem.pdac_subtyping as pds
import h5py
import numpy as np
import pandas as pd
import plotnine as gg
import scanpy as sc
from chula_stem.plotting import plot_classifier_curve
from chula_stem.r_utils import read_geo, tximport
from chula_stem.sc_rnaseq import annotate_adata_vars
from pyhere import here
from scipy import sparse
from sklearn.metrics import RocCurveDisplay, roc_curve

# %%

tx2gene = pd.read_csv(here("analyses", "data", "tx2gene.tsv"), sep="\t")
workdir = here("analyses", "pdac_subtyping")
pdata = here("analyses", "data_all", "public_data")

manifest = pd.read_csv(here(workdir, "manifest.csv"))

adata_path: Path = here(workdir, "pdac_cohort.h5ad")

if not adata_path.exists():
    adata = pds.read_manifest(manifest, tx2gene)
    adata = annotate_adata_vars(
        adata,
        savepath=here("analyses", "data", "ensembl_gene_data.csv"),
        new_index_col="hgnc_symbol",
        new_index_name="gene",
    )
    adata.X = sparse.csc_array(adata.X)
    adata.obs.loc[:, "sample"] = adata.obs_names
    adata.write_h5ad(adata_path)
else:
    adata = ad.read_h5ad(adata_path)

model = pds.MoffittBasal(threshold=0.44)  # From running last time
model.fit(adata=adata)
model.predict(inplace=True)

adata.layers["log_abundance"] = adata.layers["abundance"].copy()
sc.pp.log1p(adata, layer="log_abundance")

has_expr = adata[adata.obs["rnaseq_available"], :]
sc.pl.heatmap(
    has_expr,
    var_names=["GATA6", "SMAD4", "CDKN2A", "TP53"],
    groupby="moffitt_is_basal",
    layer="log_abundance",
    swap_axes=True,
)


def test_with_moffitt_original():
    moffitt_data = here(pdata, "GEO", "GSE71729", "GSE71729_series_matrix.txt.gz")
    moffitt = read_geo(moffitt_data)
    model = pds.MoffittBasal()
    moffitt.obs.rename(
        {"tumor_subtype_0na_1classical_2basal:ch2": "subtype"}, inplace=True, axis=1
    )
    moffitt.obs.loc[:, "is_basal-like"] = moffitt.obs["subtype"].replace(
        {"2": 1, "1": 0, "0": np.nan}
    )
    moffitt = moffitt[~moffitt.obs["is_basal-like"].isna(), :]
    model.fit(moffitt)
    model.tune(y_true="is_basal-like")
    print(f"New threshold: {model.threshold}")

    model.predict(inplace=True)
    bplot = (
        gg.ggplot(moffitt.obs, gg.aes(x="is_basal-like", y="moffitt_pr_basal"))
        + gg.geom_boxplot()
    )

    plot = plot_classifier_curve(
        curve="roc",
        y_true=moffitt.obs["is_basal-like"],
        y_score=moffitt.obs["moffitt_pr_basal"],
    )
    plot.show()
    # [2026-03-06 Fri] From running this, you got the best threshold to be 0.44
