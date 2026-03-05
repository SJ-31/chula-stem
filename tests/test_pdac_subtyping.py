#!/usr/bin/env ipython
from pathlib import Path

import anndata as ad
import chula_stem.pdac_subtyping as pds
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from chula_stem.r_utils import tximport
from chula_stem.sc_rnaseq import annotate_adata_vars
from pyhere import here
from scipy import sparse

tx2gene = pd.read_csv(here("analyses", "data", "tx2gene.tsv"), sep="\t")
workdir = here("analyses", "pdac_subtyping")

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

pds.moffitt_score(adata)
adata.obs.loc[:, "moffitt_class"] = adata.obs["moffitt_pr_basal"].apply(
    lambda x: "basal-like" if x > 0.5 else "other"
)


has_expr = adata[adata.obs["rnaseq_available"], :]

sc.pl.heatmap(
    has_expr,
    var_names=["GATA6", "SMAD4", "CDKN2A"],
    groupby="moffitt_class",
    layer="abundance",
)
