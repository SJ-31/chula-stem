#!/usr/bin/env ipython
#

import anndata as ad
import anndata2ri
import celltypist
import numpy as np
import pandas as pd
import polars as pl
import polars.selectors as cs
import scanpy as sc
from celltypist import models
from chula_stem import utils as ut
from chula_stem.sc_rnaseq import (
    annotate_marker,
    cell_assign_wrapper,
    make_marker_df,
    pca_to_leiden,
    scib_normalize,
)
from pyhere import here
from scipy import sparse

outdir = here("analyses", "output", "corneal")


def ct_annotate(
    adata: ad.AnnData, adata_final: ad.AnnData, name: str, model_file: str
) -> None:
    """Helper for annotating with celltypist"""
    to_celltypist = adata.copy()
    sc.pp.normalize_per_cell(
        to_celltypist, counts_per_cell_after=10**4
    )  # per cell typist instructions
    sc.pp.log1p(to_celltypist)
    to_celltypist.X = to_celltypist.X.toarray()
    li_model = models.Model.load(str(here("analyses", "data", model_file)))
    predictions = celltypist.annotate(
        to_celltypist, model=li_model, majority_voting=True
    ).to_adata()

    adata_final.obs["cell_type"] = predictions.obs.loc[
        adata_final.obs.index, "majority_voting"
    ]
    fig = sc.pl.umap(adata_final, color=["cell_type", "leiden"], return_fig=True)
    fig.savefig(here(outdir, f"{name}_celltypes.png"), bbox_inches="tight", dpi=500)


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
pca_to_leiden(fibro_nd, leiden_pars={"resolution": 1})

# TODO: how to quantitatively evaluate differences between these clusters in python
fibro_final, fmarkers = get_clusters(fibro_nd)
scib_normalize(  # Scran normalization on leiden clusters then log1p transform
    fibro_final, clusters=fibro_final.obs.leiden, precluster=False, cluster_method=""
)
get_hvs(fibro_final)

print("Fibro after")
print(fibro_final)

# ** Cell annotation
ct_annotate(
    fibro,
    fibro_final,
    "fibroblast",
    "Adult_Human_Skin.pkl",
)
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
pca_to_leiden(corneal, leiden_pars={"resolution": 1})
corneal_final, cmarkers = get_clusters(corneal)
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

# ** Cell annotation
markers = pl.read_csv(here("analyses", "data", "CellMarker2_human.csv"))
kept_celltypes = [
    "Fibroblast-like cell",
    "Red blood cell (erythrocyte)",
    "Blood cell",
    "White blood cell",
    "Immune cell",
]
markers = markers.filter(
    (
        (
            (pl.col("cell_type") == "Normal cell")
            & (pl.col("cell_name").is_in(kept_celltypes))
        )
        | (pl.col("tissue_class") == "Eye")
    )
    & ((pl.col("Symbol").is_not_null()) & pl.col("cell_name").is_not_null())
)
types_genes = (
    markers.select(["cell_name", "Symbol"]).group_by("cell_name").agg(pl.col("Symbol"))
)
marker_mat = make_marker_df({k: v for k, v in types_genes.iter_rows()})
cell_assign_model = here(outdir, "corneal_cellassign")
corneal_pred = cell_assign_wrapper(corneal_final, marker_mat, cell_assign_model)
fig = sc.pl.umap(
    corneal_final, color=["cell_type", "leiden"], add_outline=True, return_fig=True
)
fig.savefig(here(outdir, "corneal_celltypes.png"), bbox_inches="tight", dpi=500)

# * Get colon
# <2025-02-06 Thu> QC seems fine, no cells being discarded weirdly
# <2025-02-11 Tue> Very little LGR5 expression
colon_file = here(outdir, "colon.h5ad")
colon: ad.AnnData = ad.read_h5ad(colon_file)
print("Colon before")
print(colon)
colon = colon[(~colon.obs.discard) | (colon.obs["scDblFinder.class"] != "doublet"), :]
pca_to_leiden(colon, leiden_pars={"resolution": 1})
colon_final, clmarkers = get_clusters(colon.copy())
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
ct_annotate(
    colon, colon_final, "colon", "9_healthy_reference_AP_large_intestine_finalmodel.pkl"
)

# * Integrate

adata: ad.AnnData = ad.concat(
    [colon_final, fibro_final, corneal_final], merge="same", join="inner"
)
# Make a visual too
unintegrated = adata.copy()
pca_to_leiden(unintegrated, leiden_pars={"resolution": 1})
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
# pca_to_leiden(adata, leiden_pars={"resolution": 1})

# ** scANVI integration
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


# #  --- CODE BLOCK ---


# #  --- CODE BLOCK ---
