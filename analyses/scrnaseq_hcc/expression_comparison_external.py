#!/usr/bin/env python3

from collections.abc import Callable
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# import cellxgene_census   # BUG: can't use this on cuaim
import plotnine as gg
import scanpy as sc
import yte
from chula_stem.r_utils import edgeR_wrapper
from chula_stem.sc_rnaseq import (
    annotate_adata_vars,
    cell_assign_wrapper,
    distance_by_mads,
    mads_qc_plot_batch,
    obsm_filters,
    sc_distribution_plot,
)
from chula_stem.utils import read_existing
from loguru import logger
from matplotlib.pyplot import axis
from pyhere import here

data = here("analyses", "data_all")
outdir = here(data, "output", "HCC", "SCRNASEQ", "cohort", "compare_external")
(outdir / "qc").mkdir(exist_ok=True)
geo_dir = here(data, "public_data", "GEO")
gene_reference: Path = here("analyses/data/ensembl_gene_data.csv")

# /home/shannc/Bio_SDD/stem_synology/chula_mount/shannc/output/HCC/SCRNASEQ/cohort/compare_external

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
scvi_model_dir = outdir / "scvi_model"
max_mito_pct = 20
filter_cells_kws = {"min_genes": 200}


# CD147, GPC3 as well as all of the other genes plotted previously
# Compare with normal hepatocytes (ideally in hepatic organoids) transcriptomes.


def get_markers():  # Use only cellmarker 3.0 markers from experimental
    # evidence
    marker_df = pd.read_csv(here("analyses", "data", "human_cell_marker.txt"), sep="\t")
    marker_df = marker_df.loc[~marker_df["symbol"].isna(), :]
    marker_df = marker_df.loc[marker_df["marker_source"] == "Experiment", :]
    marker_df = marker_df.loc[marker_df["tissue_type"] == "Liver", :]
    marker_df = marker_df.loc[marker_df["tissue_class"] == "Liver", :]
    marker_df = marker_df.loc[marker_df["disease"] == "Normal", :]
    counts = marker_df["cell_name"].value_counts()
    marker_df = marker_df.loc[marker_df["cell_name"].isin(counts[counts > 5].index), :]
    return (
        marker_df.loc[:, ["cell_name", "symbol"]]
        .rename({"cell_name": "cell", "symbol": "gene"}, axis=1)
        .assign(value=1)
        .pivot_table(index="gene", columns="cell", values="value", fill_value=0)
        .astype(int)
    )


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
    for metric in ["n_genes_by_counts", "total_counts", "pct_counts_mito"]:
        qc_axes = sc.pl.violin(
            adata,
            metric,
            jitter=0.4,
            show=False,
            rotation=90.0,
            groupby="source",
        )
        qc_axes.figure.savefig(
            outdir / "qc" / f"geo_{metric}.pdf", bbox_inches="tight", dpi=500
        )
    distance_by_mads(
        adata,
        keys=["n_genes_by_counts", "total_counts", "pct_counts_mito"],
        group_keys="source",
        inplace=True,
    )
    passed_index = []
    for source in adata.obs["source"].unique():
        qc, compare = mads_qc_plot_batch(
            adata, source, batch_col="source", sample_col="sample"
        )
        (qc / compare).save(outdir / "qc" / f"{source}_mads.pdf")
        metrics_plot = (
            sc_distribution_plot(
                adata[adata.obs["source"] == source, :],
                keys=["n_genes_by_counts", "total_counts", "pct_counts_mito"],
                groupby="sample",
                fill="predicted_doublet",
                with_points=False,
            )
            + gg.ggtitle(f"Metrics for source {source}")
            + gg.theme(axis_text_x=gg.element_text(rotation=90))
        )
        metrics_plot.save(outdir / "qc" / f"{source}_main.pdf", dpi=500)
        if filter_expr := filters.get(source):
            passed_index.extend(
                list(filter_expr(adata[adata.obs["source"] == source, :]).obs_names)
            )

    adata = adata[passed_index, :]
    logger.info("{} before filtering: {}", name, adata.shape)
    adata = adata[adata.obs["pct_counts_mito"] < max_mito_pct, :]
    logger.info("{} after filtering mito: {}", name, adata.shape)
    sc.pp.filter_cells(adata, **filter_cells_kws)
    logger.info("{} after filtering cells: {}", name, adata.shape)
    return adata


# ** Organoid resources
#
# *** From GEO


def make_mads_filter_fn(spec):
    def f(adata: ad.AnnData) -> ad.AnnData:
        passed, failed = obsm_filters(adata, spec, "mads", reduction="any")
        return passed

    return f


geo_dirs = {
    "GSE166589": ["GSM5075991", "GSM5075992", "GSM5075993"],
    # GSE166589 clearly has two cell populations
    "GSE154883": "adata.h5ad",
    "GSE264261": "adata.h5ad",
    "GSE182604": "adata.h5ad",
    "GSE188541": None,
    "GSE141183": None,
    "GSE130073": None,
    "GSE207889": ["GSM6822571"],
    "GSE210059": [
        # "GSM6415980", # Excluded due to low cell count (see QC plots)
        "GSM6415982"
    ],
    # Run "expression_comparison_external.R" to get the adata files if they don't exist
    # "GSM4633164": None, # Excluded due to abnormal gene expression
}
filters: dict[
    str, Callable[[ad.AnnData], ad.AnnData]
] = {  # Filter functions to apply each data source
    "GSE130073": make_mads_filter_fn(
        {"pct_counts_mito": [None, 2], "n_genes_by_counts": [-3, None]}
    ),
    "GSE154883": make_mads_filter_fn(
        {"n_genes_by_counts": [-2, 4], "total_counts": [None, 6]}
    ),
    "GSE166589": make_mads_filter_fn({"pct_counts_mito": [None, 5]}),
    "GSE182604": make_mads_filter_fn(
        {"pct_counts_mito": [None, 3], "total_counts": [-2, 5]}
    ),
    "GSE188541": make_mads_filter_fn({"n_genes_by_counts": [-1, None]}),
    "GSE207889": make_mads_filter_fn(
        {"n_genes_by_counts": [-1.5, 3], "total_counts": [-1, 5]}
    ),
    "GSE264261": make_mads_filter_fn({"pct_counts_mito": [None, 4]}),
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
        elif isinstance(spec, str):
            adata = ad.read_h5ad(series_dir / spec)
            adata.obs["source"] = series
            if "sample" not in adata.obs.columns:
                adata.obs["sample"] = series
        else:
            tmp = []
            for sample in spec:
                sample_dir = series_dir / sample
                cur = sc.read_10x_mtx(sample_dir)
                cur.obs["sample"] = sample
                cur.obs["source"] = series
                cur.obs_names_make_unique()
                cur.var_names_make_unique()
                tmp.append(cur)
            adata = ad.concat(tmp, axis="obs", merge="first", join="outer")
        sc.pp.scrublet(adata, copy=False)
        adata.obs_names_make_unique()
        adata.var_names_make_unique()
        adatas.append(adata)
    return adatas


def get_geo(f):
    geo_adatas = collect_from_geo()
    adata: ad.AnnData = ad.concat(geo_adatas, axis="obs", merge="first", join="outer")
    adata.obs_names_make_unique()
    adata.obs["type"] = "normal"
    adata.obs["batch"] = adata.obs["source"]
    adata = annotate_filter_routine(adata, "GEO")
    adata.obs["predicted_doublet"] = adata.obs["predicted_doublet"].astype(bool)
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
    chula_doublets = []
    for p in chula.obs["patient"].unique():
        current = chula[chula.obs["patient"] == p, :]
        doublets: ad.Anndata = sc.pp.scrublet(current, copy=True)
        chula_doublets.append(
            doublets.obs.loc[:, ["doublet_score", "predicted_doublet"]]
        )
    chula.obs["source"] = "Chula"
    chula.obs["batch"] = chula.obs["flowcell"]
    chula.obs = chula.obs.loc[:, kept_cols]
    chula.obs = chula.obs.merge(
        pd.concat(chula_doublets), left_index=True, right_index=True
    )
    chula.obs["predicted_doublet"] = chula.obs["predicted_doublet"].astype(bool)
    external_geo = read_existing(external_files["geo"], get_geo, ad.read_h5ad)
    shared_cols = set(external_geo.var.columns) & set(chula.var.columns)
    chula.var = chula.var.loc[:, list(shared_cols)]
    combined = ad.concat([chula, external_geo], axis="obs", merge="first", join="outer")
    combined.var["mito"] = combined.var["mito"].astype(bool)
    sc.pp.pca(combined)
    sc.pp.neighbors(combined)
    sc.tl.umap(combined)
    sc.pp.highly_variable_genes(combined, n_top_genes=5000, flavor="seurat_v3")
    lib_size = combined.X.sum(1)
    combined.obs["size_factor"] = lib_size / np.mean(lib_size)
    cellassign_results = cell_assign_wrapper(
        combined,
        cell_markers=get_markers(),
        model_path=outdir / "cellassign",
        size_factor_key="size_factor",
    )
    combined.obs = combined.obs.merge(
        cellassign_results["pred"]["PREDICTION"], left_index=True, right_index=True
    ).rename({"PREDICTION": "cellassign_prediction"})
    if not (outdir / "cellassign").exists():
        cellassign_results["model"].save(outdir / "cellassign")
    combined = combined[
        :, combined.var["highly_variable"] | combined.var.index.isin(wanted_genes)
    ]
    combined.layers["x_norm"] = sc.pp.normalize_total(
        combined, target_sum=1e6, inplace=False
    )["X"]
    sc.pp.log1p(combined, layer="x_norm")
    combined.write_h5ad(f)
    return combined


combined: ad.AnnData = read_existing(combined_file, combine_all, ad.read_h5ad)
tmp_plots = sc.pl.umap(combined, color=["batch", "source", "type"], return_fig=True)
tmp_plots.figure.savefig(outdir / "umap_before_integration.pdf")

dotplot_check = sc.pl.dotplot(
    combined,
    groupby="sample",
    layer="x_norm",
    return_fig=True,
    var_names={
        "main": (
            "BSG",  # CD147
            "LIN28A",
            "LIN28B",
            "GPC3",  # glypican 3
            "AFP",
        ),  # Alpha-fetoprotein
        "edgeR_DE": (
            # Genes DE with edgeR
            "SDF2L1",
            "NUDT14",
            "SCAND1",
            "RPL28",
            "C4orf48",
            "HCFC1R1",
            "ABHD17A",
            "GRINA",
        ),
        "scVI_DE": ("GAL", "TAGLN", "MYL9", "RFLNB"),
    },
)
dotplot_check.savefig(outdir / "dotplot_scanpy.pdf")

# * DE analysis

# ** With scvi


def train_scvi(f):
    import scvi

    combined: ad.AnnData = read_existing(combined_file, combine_all, ad.read_h5ad)

    scvi.model.SCVI.setup_anndata(combined, batch_key="batch")
    model = scvi.model.SCVI(combined, gene_likelihood="nb")
    model.train(
        check_val_every_n_epoch=1,
        max_epochs=400,
        early_stopping=True,
        early_stopping_patience=20,
        early_stopping_monitor="elbo_validation",
    )
    model.save(f, save_anndata=False, overwrite=True)
    return model


def edgeR_de(f):
    agg_obs = combined.obs.groupby("sample").agg("first")
    adata = sc.get.aggregate(combined, by="sample", func="sum")
    adata.obs = agg_obs
    adata.X = adata.layers["sum"]
    num_de, result = edgeR_wrapper(adata, group="type", treat=False)
    result.to_csv(f, index=False)
    return result


def scvi_de(f):
    import scvi

    scvi_model: scvi.model.SCVI = read_existing(
        scvi_model_dir,
        train_scvi,
        lambda x: scvi.model.SCVI.load(x, adata=combined),
    )
    result: pd.DataFrame = scvi_model.differential_expression(
        groupby="sample",
        mode="change",
        batch_correction=True,
        idx1=combined.obs["type"] == "tumor",
        idx2=combined.obs["type"] == "normal",
    ).reset_index(names="gene")
    result.to_csv(f, index=False)
    return result


edger_de_results = read_existing(outdir / "edger_de_results.csv", edgeR_de, pd.read_csv)
scvi_de_results = read_existing(outdir / "scvi_de_results.csv", scvi_de, pd.read_csv)
