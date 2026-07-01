#!/usr/bin/env python3

from pathlib import Path

import anndata as ad

# import cellxgene_census   # BUG: can't use this on cuaim
import plotnine as gg
import scanpy as sc
import yte
from chula_stem.sc_rnaseq import annotate_adata_vars, distance_by_mads
from chula_stem.utils import read_existing
from loguru import logger
from pyhere import here

data = here("analyses", "data_all")
outdir = here(data, "output", "HCC", "SCRNASEQ", "cohort", "compare_external")
geo_dir = here(data, "public_data", "GEO")
gene_reference: Path = here("analyses/data/ensembl_gene_data.csv")

wanted_genes = [
    "BSG",  # CD147
    "LIN28A",
    "LIN28B",
    "GPC3",  # glypican 3
    "AFP",  # Alpha-fetoprotein
]

external_files: dict[str, Path] = {
    "cellxgene": outdir / "external_cellxgene.h5ad",
    "geo": outdir / "external_geo.h5ad",
}
combined_file = outdir / "combined.h5ad"
max_mito_pct = 20
filter_cells_kws = {"min_genes": 200}

# CD147, GPC3 as well as all of the other genes plotted previously
# Compare with normal hepatocytes (ideally in hepatic organoids) transcriptomes.


# * Data import


def annotate_filter_routine(adata: ad.AnnData, name: str) -> ad.AnnData:
    adata = annotate_adata_vars(
        adata,
        merge_on="gene_ids",
        savepath=gene_reference,
        meta_id_col="ensembl_gene_id",
    )
    logger.info("Succesfully annotated {} adata\n{}", name, adata.var)
    sc.pp.calculate_qc_metrics(adata, inplace=True, qc_vars=["mito"])
    logger.info("{} before filtering: {}", name, adata.shape)
    adata = adata[adata.obs["pct_counts_mito"] < max_mito_pct, :]
    logger.info("{} after filtering mito: {}", name, adata.shape)
    sc.pp.filter_cells(adata, **filter_cells_kws)
    logger.info("{} after filtering cells: {}", name, adata.shape)
    return adata


# ** Organoid resources
#
# *** From GEO

geo_dirs = {
    "GSE166589": ["GSM5075991", "GSM5075992", "GSM5075993"],
    "GSE154883": "adata.h5ad",
    # Run "expression_comparison_external.R" to get the file above
    "GSM4633164": None,
}


def collect_from_geo() -> list[ad.AnnData]:
    adatas = []
    for series, spec in geo_dirs.items():
        series_dir = geo_dir / series
        if spec is None:
            adata = sc.read_10x_mtx(series_dir)
            adata.obs["source"] = series
            if "sample" not in adata.obs.columns:
                adata.obs["sample"] = series
            adatas.append(adata)
        elif isinstance(spec, str):
            adata = ad.read_h5ad(series_dir / spec)
            adata.obs["source"] = series
            if "sample" not in adata.obs.columns:
                adata.obs["sample"] = series
            adatas.append(adata)
        else:
            for sample in spec:
                sample_dir = series_dir / sample
                adata = sc.read_10x_mtx(sample_dir)
                adata.obs["sample"] = sample
                adata.obs["source"] = series
                adatas.append(adata)
    return adatas


def get_geo(f):
    geo_adatas = collect_from_geo()
    adata: ad.AnnData = ad.concat(geo_adatas, axis="obs", merge="first", join="outer")
    adata.obs["type"] = "normal"
    adata.obs["batch"] = adata.obs["source"]
    adata.obs_names_make_unique()
    adata = annotate_filter_routine(adata, "GEO")
    adata.write_h5ad(f)
    return adata


# *** From cellxgene

# [2026-06-30 Tue] No organoid data with hepatocytes available on cellxgene
# closest thing are hepatoblasts

# if not external_files["cellxgene"].exists():
#     with cellxgene_census.open_soma() as census:
#         normal_tissues = cellxgene_census.get_anndata(
#             census=census,
#             organism="Homo sapiens",
#             census_version="2025-11-08",
#             obs_value_filter="tissue == 'hepatocyte'",
#             obs_column_names={
#                 "obs": ["assay", "cell_type", "tissue", "disease", "tissue_general"]
#             },
#         )
# else:
#     external_cellxgene = ad.read_h5ad(external_files["cellxgene"], backed=True)


# ** Primary tissue
# Use cellxgene to get hepatocytes from normal tissues

# with cellxgene_census.open_soma() as census:
#     external_tissue = cellxgene_census.get_anndata(
#         census=census,
#         organism="Homo sapiens",
#         census_version="2025-11-08",
#         obs_value_filter="cell_type == 'hepatocyte'",
#         obs_column_names={
#             "obs": [
#                 "cell_type",
#                 "tissue",
#                 "tissue_general",
#                 "dataset_id",
#             ]
#         },
#     )

# * Combine with all


def combine_all(f):
    kept_cols = [
        "sample",
        "source",
        "type",
        "batch",
        "log1p_n_genes_by_counts",
        "total_counts",
        "log1p_total_counts",
        "pct_counts_in_top_50_genes",
        "pct_counts_in_top_100_genes",
        "pct_counts_in_top_200_genes",
        "pct_counts_in_top_500_genes",
        "total_counts_mito",
        "log1p_total_counts_mito",
        "pct_counts_mito",
        "lshift_size_factors",
    ]
    chula = ad.read_h5ad(here(data, "output", "HCC", "SCRNASEQ/cohort/filtered.h5ad"))
    chula = chula[
        (chula.obs["treatment"] == "control") & (chula.obs["type"] == "tumor"), :
    ]
    chula.obs["source"] = "Chula"
    chula.obs["batch"] = chula.obs["flowcell"]
    chula.obs = chula.obs.loc[:, kept_cols]
    external_geo = read_existing(external_files["geo"], get_geo, ad.read_h5ad)
    shared_cols = set(external_geo.var.columns) & set(chula.var.columns)
    chula.var = chula.var.loc[:, list(shared_cols)]
    combined = ad.concat([chula, external_geo], axis="obs", merge="first", join="outer")
    combined.var["mito"] = combined.var["mito"].astype(bool)
    sc.pp.pca(combined)
    sc.pp.neighbors(combined)
    sc.tl.umap(combined)
    combined.layers["x_norm"] = sc.pp.normalize_total(
        combined, target_sum=1e6, inplace=False
    )["X"]
    sc.pp.log1p(combined, layer="x_norm")
    combined.write_h5ad(f)
    return combined


combined = read_existing(
    combined_file, combine_all, lambda x: ad.read_h5ad(x, backed=True)
)
tmp_plots = sc.pl.umap(combined, color=["batch", "source", "type"], return_fig=True)
tmp_plots.figure.savefig(outdir / "umap_before_integration.pdf")
