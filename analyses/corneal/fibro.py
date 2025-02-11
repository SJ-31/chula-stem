#!/usr/bin/env ipython
#

import json

import anndata as ad
import anndata2ri
import celltypist
import numpy as np
import pandas as pd
import polars as pl
import polars.selectors as cs
import scanpy as sc
import scanpy.external as sce
import scvi
from celltypist import models
from chula_stem import utils as ut
from chula_stem.sc_rnaseq import cell_assign_wrapper, make_marker_df
from pyhere import here
from scipy import sparse
from scvi.external import CellAssign

outdir = here("analyses", "output", "corneal")


# Modified version of normalize from scib that won't throw error with zeros
def scib_normalize(
    adata,
    min_mean=0.1,
    log=True,
    precluster=True,
    clusters=(),
    cluster_method="louvain",
    sparsify=True,
):
    import rpy2.rinterface_lib.callbacks
    import rpy2.rinterface_lib.embedded
    import rpy2.robjects as ro

    # Check for 0 count cells
    if np.any(adata.X.sum(axis=1) == 0):
        print("WARNING: found 0 count cells in the AnnData object")
        print("Filtering cells...")
        sc.pp.filter_cells(adata, min_counts=10, inplace=True)

    # Check for 0 count genes
    if np.any(adata.X.sum(axis=0) == 0):
        sc.pp.filter_genes(adata, min_counts=1, inplace=True)
        print("WARNING: found 0 count genes in the AnnData object")

    if sparsify:
        # massive speedup when working with sparse matrix
        if not sparse.issparse(adata.X):  # quick fix: HVG doesn't work on dense matrix
            adata.X = sparse.csr_matrix(adata.X)
    try:
        ro.r("library(scran)")
    except rpy2.rinterface_lib.embedded.RRuntimeError as ex:
        raise ImportError("scran library not installed")

    anndata2ri.activate()

    # keep raw counts
    adata.layers["counts"] = adata.X.copy()

    is_sparse = False
    x = adata.X.T
    # convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
    if sparse.issparse(x):
        is_sparse = True

        if x.nnz > np.power(2, 31) - 1:
            x = x.tocoo()
        else:
            x = x.tocsc()

    ro.globalenv["data_mat"] = x
    if len(clusters) > 0:
        precluster = False
    if precluster:
        # Preliminary clustering for differentiated normalisation
        adata_pp = adata.copy()
        sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6, min_counts=1)
        sc.pp.log1p(adata_pp)
        sc.pp.pca(adata_pp, n_comps=15, svd_solver="arpack")
        sc.pp.neighbors(adata_pp)
        if cluster_method == "louvain":
            sc.tl.louvain(adata_pp, key_added="groups", resolution=0.5)
        elif cluster_method == "leiden":
            sc.tl.leiden(adata_pp, key_added="groups", resolution=0.5)
        else:
            raise NotImplementedError(
                "Choose `cluster_method` from 'louvain', 'leiden'"
            )

        ro.globalenv["input_groups"] = adata_pp.obs["groups"]
        size_factors = ro.r(
            "sizeFactors("
            "   computeSumFactors("
            "       SingleCellExperiment(list(counts=data_mat)),"
            "       clusters = input_groups,"
            f"       min.mean = {min_mean}"
            "   )"
            ")"
        )

        del adata_pp
    elif len(clusters) > 0:
        ro.globalenv["clusters"] = clusters
        size_factors = ro.r(
            f"""
            sizeFactors(computeSumFactors(SingleCellExperiment(
            list(counts=data_mat)), min.mean = {min_mean}, clusters=clusters))
            """
        )
    else:
        size_factors = ro.r(
            f"""
            sizeFactors(computeSumFactors(SingleCellExperiment(
            list(counts=data_mat)), min.mean = {min_mean}))
            """
        )

    # modify adata
    adata.obs["size_factors"] = size_factors
    adata.X /= adata.obs["size_factors"].values[:, None]
    if log:
        print("Note! Performing log1p-transformation after normalization.")
        sc.pp.log1p(adata)
    else:
        print("No log-transformation performed after normalization.")

    if is_sparse:
        # convert to sparse, bc operation always converts to dense
        adata.X = sparse.csr_matrix(adata.X)

    adata.raw = adata  # Store the full data set in 'raw' as log-normalised data for statistical testing

    # Free memory in R
    ro.r("rm(list=ls())")
    ro.r("gc()")

    anndata2ri.deactivate()


def pca_cluster(
    adata, pca_pars=None, neighbor_pars=None, umap_pars=None, leiden_pars=None
):
    if "X_pca" not in adata.obsm:
        ut.do_call(lambda **x: sc.pp.pca(adata, **x), pca_pars)
    ut.do_call(lambda **x: sc.pp.neighbors(adata, **x), neighbor_pars)
    ut.do_call(lambda **x: sc.tl.umap(adata, **x), umap_pars)
    ut.do_call(lambda **x: sc.tl.leiden(adata, **x), leiden_pars)


def annotate_marker(
    adata: ad.AnnData, marker_genes: list, gene_col: str = None
) -> pd.DataFrame:
    expr: pd.DataFrame = (
        sc.get.obs_df(adata, keys=marker_genes, gene_symbols=gene_col) > 0
    )
    for c in expr.columns:
        adata.obs[f"has_{c}"] = expr[c]
    return expr.agg("sum")


def get_clusters(adata) -> tuple[ad.AnnData, pl.DataFrame]:
    """Filter adata object to contain desired clusters
    <2025-02-10 Mon> will keep clusters that have >= 40% cells with nonzero LGR5 expression
        Confirm with P'Ta if this is what she wants, because you need to use clusters when
        doing pseudobulk DGE
    """
    _ = annotate_marker(adata, MARKERS)
    marker_df = pl.from_pandas(adata.obs).select(cs.starts_with("has_"))
    grouped: pd.DataFrame = adata.obs.groupby("leiden").agg(
        has_LGR5=("has_LGR5", "sum"), size=("has_LGR5", "size")
    )
    print(grouped)
    threshold = 0.4
    pass_threshold = grouped.query("has_LGR5/size >= @threshold").index.to_list()
    passed = adata[adata.obs.leiden.isin(pass_threshold), :]
    return adata, marker_df  # <2025-02-11 Tue> temporarily just return the whole object


# while we figure out what to filter on


MARKERS = ["LGR5", "OR4F5", "TET1", "TET2", "TET3"]


def get_hvs(adata) -> None:
    sc.pp.highly_variable_genes(
        adata, flavor="cell_ranger", n_top_genes=2000, subset=True
    )


# * Get fibro
# <2025-02-06 Thu> Cells with high mitochondrial counts are unusually abundant
# Many of these also have high total counts and detected genes, so unlikely to be
# damaged cells
# Will change "discard" column to reflect this: do not discard if BOTH sum and detection
# filters pass

fibro_file = here(outdir, "fibro.h5ad")
fibro: ad.AnnData = ad.read_h5ad(fibro_file)
print("Fibro before")
print(fibro)
# #  --- CODE BLOCK ---
# Discard low-quality and doublets
#
fibro.obs.discard = ~fibro.obs.discard_reason.isin(
    ["kept", "high_subsets_mito_percent"]
)

fibro = fibro[~fibro.obs.discard, :]
fibro_nd = fibro[fibro.obs["scDblFinder.class"] != "doublet", :]
# pca_cluster(fibro, leiden_pars={"resolution": 1})
pca_cluster(fibro_nd, leiden_pars={"resolution": 1})

fig = sc.pl.umap(fibro_nd, color="leiden", add_outline=True, return_fig=True)
fig.savefig(here(outdir, "fibro_leiden.png"), bbox_inches="tight", dpi=500)

# TODO: how to quantitatively evaluate differences between these clusters in python
fibro_final, fmarkers = get_clusters(fibro_nd)
scib_normalize(  # Scran normalization on leiden clusters then log1p transform
    fibro_final, clusters=fibro_final.obs.leiden, precluster=False, cluster_method=""
)
get_hvs(fibro_final)

print("Fibro after")
print(fibro_final)

# ** Cell annotation
# <2025-02-11 Tue> Based on the marker genes provided by original authors, using

with open(here("analyses", "corneal", "fibro_types.json"), "r") as f:
    fibro_markers = make_marker_df(json.load(f))


fibro_pred = cell_assign_wrapper(fibro_final, fibro_markers)

sc.pl.umap(fibro_final, color=["cell_type", "leiden"], ncols=1)

# #  --- CODE BLOCK ---

# * Get corneal
# <2025-02-07 Fri> Same deal as fibro sample: cells with high mitochondrial counts but
# high total counts and detected genes are present
#
# <2025-02-11 Tue> This has no LGR5 expression...
corneal_file = here(outdir, "corneal.h5ad")
corneal: ad.AnnData = ad.read_h5ad(corneal_file)
print("Corneal before")
print(corneal)
corneal.obs.discard = ~corneal.obs.discard_reason.isin(
    ["kept", "high_subsets_mito_percent"]
)
corneal = corneal[~corneal.obs.discard, :]
corneal = corneal[corneal.obs["scDblFinder.class"] != "doublet", :]
pca_cluster(corneal, leiden_pars={"resolution": 1})
corneal_final, cmarkers = get_clusters(corneal)
fig = sc.pl.umap(corneal, color="leiden", add_outline=True, return_fig=True)
fig.savefig(here(outdir, "corneal_leiden.png"), bbox_inches="tight", dpi=500)
scib_normalize(
    corneal_final,
    clusters=corneal_final.obs.leiden,
    precluster=False,
    cluster_method="",
)
get_hvs(corneal_final)
print("Corneal after")
print(corneal_final)
# #  --- CODE BLOCK ---

# * Get colon
# <2025-02-06 Thu> QC seems fine, no cells being discarded weirdly
# <2025-02-11 Tue> Very little LGR5 expression
colon_file = here(outdir, "colon.h5ad")
colon: ad.AnnData = ad.read_h5ad(colon_file)
print("Colon before")
print(colon)
colon = colon[(~colon.obs.discard) | (colon.obs["scDblFinder.class"] != "doublet"), :]
pca_cluster(colon, leiden_pars={"resolution": 1})
colon_final, clmarkers = get_clusters(colon.copy())
fig = sc.pl.umap(colon_final, color="leiden", add_outline=True, return_fig=True)
fig.savefig(here(outdir, "colon_leiden.png"), bbox_inches="tight", dpi=500)
scib_normalize(
    colon_final,
    clusters=colon_final.obs.leiden,
    precluster=False,
    cluster_method="",
)
get_hvs(colon_final)
print("Colon after")
print(colon_final)

# ** Cell annotation
to_celltypist = colon.copy()
sc.pp.normalize_per_cell(
    to_celltypist, counts_per_cell_after=10**4
)  # per cell typist instructions
sc.pp.log1p(to_celltypist)
to_celltypist.X = to_celltypist.X.toarray()
model_file = "9_healthy_reference_AP_large_intestine_finalmodel.pkl"
li_model = models.Model.load(str(here("analyses", "data", model_file)))
predictions = celltypist.annotate(
    to_celltypist, model=li_model, majority_voting=True
).to_adata()

colon_final.obs["cell_type"] = predictions.obs.loc[
    colon_final.obs.index, "majority_voting"
]
fig = sc.pl.umap(colon_final, color=["cell_type", "leiden"], return_fig=True)
fig.savefig(here(outdir, "colon_types.png"), bbox_inches="tight", dpi=500)


# * Integrate

adata: ad.AnnData = ad.concat(
    [colon_final, fibro_final, corneal_final], merge="same", join="inner"
)
# Make a visual too
unintegrated = adata.copy()
pca_cluster(unintegrated, leiden_pars={"resolution": 1})
fig = sc.pl.umap(fibro_nd, color="source", add_outline=True, return_fig=True)
fig.savefig(here(outdir, "all_leiden.png"), bbox_inches="tight", dpi=500)
print("Integration")
print(adata)

marker_df = pl.concat([cmarkers.sum(), fmarkers.sum(), clmarkers.sum()]).with_columns(
    source=pl.Series(["Corneal", "Fibro", "Colon"]),
    n_cells=pl.Series([a.shape[0] for a in [corneal, fibro, colon]]),
)
marker_df.write_csv(here(outdir, "marker_counts.csv"))

# <2025-02-11 Tue> But batch is FULLY confounded with source. What to do?
# batch correction might remove all of the variation
# So need to use a method that handles this situation well - possibly scGen and scANVI
# sc.pp.pca(adata)
# sce.pp.scanorama_integrate(adata, "source", verbose=1)
# pca_cluster(adata, leiden_pars={"resolution": 1})

# ** scANVI integration
# Requires cell annotations... <2025-02-11 Tue> just annotate from what the original authors did
#
# adata.raw = adata  # keep full dimension safe
# scvi.model.SCVI.setup_anndata(adata, batch_key="source")

# model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
# model.train()

# scanvi_model = scvi.model.SCANVI.from_scvi_model(
#     model,
#     adata=adata,
#     labels_key="cell_type",
#     unlabeled_category="Unknown",
# )
# scanvi_model.train(max_epochs=20, n_samples_per_label=100)
# SCANVI_LATENT_KEY = "X_scANVI"
# adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)
# sc.pl.umap(
#     adata,
#     color=["cell_type"],
#     frameon=False,
#     ncols=1,
# )
