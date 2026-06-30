#!/usr/bin/env python3

from pathlib import Path

import anndata as ad
import cellxgene_census
import plotnine as gg
import scanpy as sc
import yte
from chula_stem.sc_rnaseq import annotate_adata_vars, distance_by_mads
from loguru import logger
from pyhere import here

data = here("analyses", "data_all")
outdir: Path = here("analyses", "output", "scrnaseq_hcc")
geo_dir = here(data, "public_data", "GEO")
gene_reference: Path = here("analyses/data/ensembl_gene_data.csv")

wanted_genes = [
    "BSG",  # CD147
    "LIN28A",
    "LIN28B",
    "GPC3",  # glypican 3
    "AFP",  # Alpha-fetoprotein
]

combined = ad.read_h5ad(
    here(data, "output", "HCC", "SCRNASEQ/cohort/filtered.h5ad"), backed=True
)
combined = combined[combined.obs["treatment"] == "control", :]

external_geo_file: Path = here(
    data, "output", "HCC", "SCRNASEQ", "cohort", "external_combined.h5ad"
)

max_mito_pct = 20
filter_cells_kws = {"min_genes": 200}


# CD147, GPC3 as well as all of the other genes plotted previously
# Compare with normal hepatocytes (ideally in hepatic organoids) transcriptomes.


# Use cellxgene to get hepatocytes from normal tissues

# * Data import

# ** Organoid resources
#
# *** From GEO

geo_dirs = {
    "GSE166589": ["GSM5075991", "GSM5075992", "GSM5075993"],
    # "GSE154883": None, can't use, is multiplexed with diseased samples
    "GSM4633164": None,
}


def collect_from_geo() -> list[ad.AnnData]:
    adatas = []
    for series, samples in geo_dirs.items():
        series_dir = geo_dir / series
        if samples is None:
            adata = sc.read_10x_mtx(series_dir)
            adata.obs["source"] = series
            adatas.append(adata)
        else:
            for sample in samples:
                sample_dir = series_dir / sample
                adata = sc.read_10x_mtx(sample_dir)
                adata.obs["source"] = sample
                adatas.append(adata)
    return adatas


if not external_geo_file.exists():
    geo_adatas = collect_from_geo()
    external_geo: ad.AnnData = ad.concat(geo_adatas, axis="obs", merge="first")
    external_geo.obs_names_make_unique()
    external_geo = annotate_adata_vars(
        external_geo,
        merge_on="gene_ids",
        savepath=gene_reference,
        meta_id_col="ensembl_gene_id",
    )
    logger.info("Succesfully annotated GEO adata\n{}", external_geo.var)
    sc.pp.calculate_qc_metrics(external_geo, inplace=True, qc_vars=["mito"])
    logger.info("GEO before filtering: {}", external_geo.shape)
    external_geo = external_geo[external_geo.obs["pct_counts_mito"] < max_mito_pct, :]
    logger.info("GEO after filtering mito: {}", external_geo.shape)
    sc.pp.filter_cells(external_geo, **filter_cells_kws)
    logger.info("GEO after filtering cells: {}", external_geo.shape)
    external_geo.write_h5ad(external_geo_file)
else:
    external_geo = ad.read_h5ad(external_geo_file, backed=True)


# *** From cellxgene

# with cellxgene_census.open_soma() as census:
#     normal_tissues = cellxgene_census.get_anndata(
#         census=census,
#         organism="Homo sapiens",
#         census_version="2025-11-08",
#         obs_value_filter="cell_type == 'hepatocyte'",
#         obs_column_names={
#             "obs": ["assay", "cell_type", "tissue", "tissue_general", "suspension_type"]
#         },
#     )
