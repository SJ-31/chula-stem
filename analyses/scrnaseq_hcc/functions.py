#!/usr/bin/env python

import re
from functools import reduce
from pathlib import Path
from typing import Literal

import anndata as ad
import marimo as mo
import numpy as np
import pandas as pd
import plotnine as gg
import polars as pl
import scanpy as sc
import yaml
from chula_stem.sc_rnaseq import annotate_adata_vars, annotate_marker, distance_by_mads
from chula_stem.utils import read_existing
from loguru import logger

# * Data prep/retrieval


def prepare_data(file, env):
    droot = Path(env["data_root"])
    tmp = []
    manifest = pl.read_csv(env["manifest"])
    target = "filtered_feature_bc_matrix.h5"
    for row in manifest.iter_rows(named=True):
        patient = row["sample_name"]
        stype = row["type"]
        suffix = f"{row['type']}-{row['treatment']}"
        data_path: Path = (
            droot / patient / "processed" / f"cellranger_{suffix}/{target}"
        )
        if not data_path.exists():
            logger.warning("The files for sample {} ({}) don't exist", patient, suffix)
            continue
        current: ad.AnnData = sc.read_10x_h5(data_path)
        current.obs_names_make_unique()
        current.var_names_make_unique()
        current.obs.loc[:, "sample"] = f"{patient}_{suffix}"
        current.obs.loc[:, "patient"] = patient
        current.obs.loc[:, "flowcell"] = row["flowcell"]
        current.obs.loc[:, "type"] = stype
        current.obs.loc[:, "treatment"] = row["treatment"]
        tmp.append(current)
    if extra := env["files"].get("extras"):
        for key, path in extra.items():
            patient, stype, treatment = key.split("|")
            if not Path(path).exists():
                logger.warning("Path to extra file {} does not exist", path)
                continue
            current = ad.read_h5ad(path)
            current.obs.loc[:, "patient"] = patient
            current.obs.loc[:, "sample"] = f"{patient}_{stype}-{treatment}"
            current.obs.loc[:, "type"] = stype
            current.obs.loc[:, "treatment"] = treatment
            tmp.append(current)
    adata = ad.concat(
        tmp, merge="first", index_unique="-", join="outer", uns_merge="same"
    )
    adata = adata[:, ~adata.var["gene_ids"].isna()]
    annotate_adata_vars(adata, "gene_ids", savepath=Path(env["gene_reference"]))
    annotate_marker(adata, marker_genes=env["markers"]["main"], gene_col="hgnc_symbol")
    sc.pp.calculate_qc_metrics(adata, inplace=True, qc_vars=["mito"])
    # TODO: should you change the grouping col?
    # review after looking at the data unintegrated
    distance_by_mads(
        adata,
        ["total_counts", "n_genes_by_counts", "pct_counts_mito"],
        group_keys=None,
        inplace=True,
    )
    if env.get("test", False):
        sc.pp.subsample(adata, fraction=0.1)
        adata = adata[:, :1000]
    with_normalized_layers(adata, env=env)
    # PCA with permissive setting for HVGs to speed things up
    hvgs: pd.DataFrame = sc.pp.highly_variable_genes(
        adata,
        inplace=False,
        layer="lshift_normalized",
        n_top_genes=5000,
        batch_key="sample",
    )
    sc.pp.pca(
        adata, n_comps=100, layer="lshift_normalized", mask_var=hvgs["highly_variable"]
    )
    adata.write_h5ad(file)
    adata = ad.read_h5ad(file, backed=True)
    return adata


def add_cfs_community(
    adata: ad.AnnData,
    old_clusters: Path | str,
    cfs_community_file: Path | str,
    column: str = "cfs_community",
):
    """
    Reassign clusters using the results of a previous call of ClusterFoldSimilarity
    """
    if isinstance(old_clusters, Path):
        cells2old_clusters = pd.read_csv(
            old_clusters, names=["cell_id", "tmp_cluster"]
        ).set_index("cell_id")
        old = cells2old_clusters[adata.obs_names]
    else:
        old = adata.obs[old_clusters]
    communities = pd.read_csv(cfs_community_file)
    lookup = {(s, g): c for s, g, c in communities.iterrows()}
    adata.obs = adata.obs.merge(old, left_index=True, right_index=True)
    reassigned = adata.obs["sample"].combine(
        adata.obs["tmp_cluster"], lambda x, y: lookup[(x, y)]
    )
    adata.obs.loc[:, column] = reassigned


def get_gs_and_cell_markers(env: dict) -> tuple[dict, dict]:
    cell_markers = env["markers"] or {}
    gene_sets = env["gene_sets"] or {}
    for update_to, file_list in zip(
        (cell_markers, gene_sets), ("marker_files", "gene_set_files")
    ):
        for file in env[file_list]:
            with open(file, "r") as f:
                update_to.update(yaml.safe_load(f))
    return gene_sets, cell_markers


def with_normalized_layers(adata: ad.AnnData, env: dict) -> None:
    "Add additional normalization layers to `adata`, retaining the raw data in X"
    cfg = env["normalization"]
    # BUG: this wasn't working
    lshift = sc.pp.normalize_total(adata, inplace=False, **cfg["normalize_total"])
    adata.layers["lshift_normalized"] = lshift["X"]
    adata.layers["lshift_size_factors"] = lshift["norm_factor"]
    sc.pp.log1p(adata, layer="lshift_normalized")
    # pn = pooled_normalization(adata, inplace=False, **cfg["pooled_normalization"])
    # adata.layers["scran_normalized"] = pn
    # TODO: think about doing it within batches
    # if cfg["normalize_total"].get("within_batch", False):


def data_import(env: dict) -> ad.AnnData:
    adata = read_existing(
        Path(env["files"]["combined"]),
        lambda x: prepare_data(x, env=env),
        lambda x: ad.read_h5ad(x, backed=True),
    )
    print(f"Samples present: {set(adata.obs["sample"])}")
    return adata


# * Plotting


def save_dotplots(
    adata: ad.AnnData, markers, grouping_keys, outdir: Path, env: dict
) -> None:
    sc.pl.DotPlot.DEFAULT_COLORMAP = "Blues"
    for key in grouping_keys:
        plot = sc.pl.DotPlot(
            adata, var_names=markers, groupby=key, **env["dotplot_kws"]
        )
        plot.fig.savefig(outdir / f"{key}.pdf")


def qc_plot_patient(adata: ad.AnnData, patient, thresholds=[1, 3, 5], line_alpha=0.5):
    adata = adata[adata.obs["patient"] == patient, :].copy()
    adata.obs.loc[:, "sample"] = adata.obs["sample"].str.replace(".*_", "")
    plots = {}
    nudge_by = len(set(adata.obs["sample"]))
    for i, col in enumerate(["total_counts", "pct_counts_mito", "n_genes_by_counts"]):
        plot = (
            gg.ggplot(
                adata.obs,
                gg.aes(x="sample", y=adata.obsm["mads"][col], fill="sample"),
            )
            + gg.geom_violin(position=gg.position_nudge(x=nudge_by))
            + gg.theme(figure_size=(15, 20), axis_title_x=gg.element_blank())
            + gg.guides(alpha="none")
            + gg.geom_point(
                gg.aes(color=col, alpha=0.1),
                position=gg.position_jitter(),
                fill=None,
                size=0.3,
            )
            + gg.scale_color_continuous("cool")
            + gg.ylab(f"{col} (MADs distance from median)")
        )
        if i == 0:
            plot = plot + gg.ggtitle(f"Patient: {patient}")
        if i != 2:
            plot = plot + gg.theme(
                axis_text_x=gg.element_blank(), axis_ticks_major_x=gg.element_blank()
            )
        col_max, col_min = adata.obsm["mads"][col].max(), adata.obsm["mads"][col].min()
        for threshold in thresholds:
            if col_max > threshold:
                plot = plot + gg.geom_hline(yintercept=threshold, alpha=line_alpha)
            if col_min < -threshold:
                plot = plot + gg.geom_hline(yintercept=-threshold, alpha=line_alpha)
        plots[i] = plot
    cplot = (
        gg.ggplot(
            pd.concat([adata.obsm["mads"], adata.obs.loc[:, ["sample"]]], axis=1),
            gg.aes(x="total_counts", y="n_genes_by_counts", color="pct_counts_mito"),
        )
        + gg.geom_point(position=gg.position_jitter())
        + gg.scale_color_continuous("coolwarm")
        + gg.facet_wrap("sample", scales="free")
        + gg.theme(figure_size=(15, 20))
        + gg.ggtitle("Relationship between QC variables (expressed as MADs distances)")
    )
    y_max, y_min = (
        adata.obsm["mads"]["n_genes_by_counts"].max(),
        adata.obsm["mads"]["n_genes_by_counts"].min(),
    )
    x_max, x_min = (
        adata.obsm["mads"]["total_counts"].max(),
        adata.obsm["mads"]["total_counts"].min(),
    )
    for threshold in thresholds:
        if y_max > threshold:
            cplot = cplot + gg.geom_hline(
                yintercept=threshold, alpha=line_alpha, linetype="dashed"
            )
        if y_min < -threshold:
            cplot = cplot + gg.geom_hline(
                yintercept=-threshold, alpha=line_alpha, linetype="dashed"
            )
        if x_max > threshold:
            cplot = cplot + gg.geom_vline(
                xintercept=threshold, alpha=line_alpha, linetype="dotted"
            )
        if x_min < -threshold:
            cplot = cplot + gg.geom_vline(
                xintercept=-threshold, alpha=line_alpha, linetype="dotted"
            )
    return plots[0] / plots[1] / plots[2], cplot


def add_saved_dr(adata: ad.AnnData, env: dict, dir_key: str = "unintegrated") -> None:
    dr_dir: Path = Path(env["DR"]["outdirs"][dir_key])
    methods = env["DR"]["methods"].keys()
    for d in dr_dir.iterdir():
        if not d.is_dir() and d.stem not in methods:
            continue
        method: str = d.stem
        for efile in d.iterdir():
            embeddings = np.load(efile)
            hp_val = efile.stem.split("_", 1)[1]
            adata.obsm[f"X_{method}_{hp_val}"] = embeddings


def plot_dr(
    adata: ad.AnnData,
    key: str,
    color_by: list[str] | str | None = None,
    prefix: str | None = None,
    axis_suffixes=("1", "2"),
    **theme_kws,
) -> list[gg.ggplot] | gg.ggplot:
    prefix = prefix or key
    x, y = f"{prefix}{axis_suffixes[0]}", f"{prefix}{axis_suffixes[1]}"
    vals = pd.DataFrame(adata.obsm[key], columns=[x, y], index=adata.obs_names)
    color_by = [color_by] if isinstance(color_by, str) else color_by
    if color_by is not None:
        vals = pd.concat([vals, adata.obs.loc[:, color_by]], axis="columns")
    plot = gg.ggplot(vals, gg.aes(x=x, y=y))
    if color_by is None:
        return plot + gg.geom_point(size=0.3)
    plots = []
    for i, col in enumerate(color_by):
        cur = plot + gg.geom_point(gg.aes(color=col), size=0.3) + gg.theme(**theme_kws)
        if i > 0:
            cur = cur + gg.theme(axis_title_y=gg.element_blank())
        plots.append(cur)
    if len(plots) == 1:
        return plots[0]
    return plots


def make_dr_slider(
    adata: ad.AnnData,
    dr_name: Literal["t-sne", "umap"],
    outdir: Path,
    color_by: list[str] | None = None,
    theme_kws: dict | None = None,
    ext: str = "pdf",
    **save_kws,
):
    """
    Return a slider object and image call
    for displaying dimensionality reduction plots in marimo,
    varying by a hyperparameter value
    """
    theme_kws = theme_kws or {}
    hp_vals = []
    val2keys, val2plots = {}, {}
    for k in adata.obsm.keys():
        if k.startswith(f"X_{dr_name}"):
            val = int(re.findall("_([0-9]+$)", k)[0])
            val2keys[val] = k
            hp_vals.append(val)

    outdir.mkdir(exist_ok=True)

    for k, v in val2keys.items():
        fname = outdir / f"{k}.{ext}"
        if fname.exists():
            continue
        plots = plot_dr(
            adata, key=v, prefix=dr_name.upper(), color_by=color_by, **theme_kws
        )
        if isinstance(plots, list):
            plots = reduce(lambda x, y: x | y, plots)
        val2plots[k] = plots
        plots.save(fname, verbose=False, **save_kws)

    slider = mo.ui.slider(steps=sorted(hp_vals), label=dr_name, show_value=True)

    def display_call(v):
        if ext == "pdf":
            return mo.pdf(outdir / f"{v}.pdf")
        return mo.image(outdir / f"{v}.{ext}")

    return slider, display_call


# * QC


def mads_filter_outliers(
    adata: ad.AnnData,
    filters: dict[str, tuple[float | None, float | None]],
    mads_key: str = "mads",
    reduction: Literal["all", "any"] = "all",
) -> tuple[ad.AnnData, ad.AnnData]:
    """Remove outliers based on mads thresholds

    Parameters
    ----------
    filters : dict
            COLNAME: [<OUTLIER_LOW>, <OUTLIER_HIGH>]
    A cell is determined to be an outlier if its distance from the median in terms of MADs falls below OUTLIER_LOW or is greater than OUTLIER_HIGH. They can be set to null not
    to do filtering in that direction

    Returns
    -------
    tuple of adata_passing, adata_failed


    Notes
    -----
    This function should be used after calling `distance_by_mads`

    """
    mads_df: pd.DataFrame = adata.obsm[mads_key]
    masks = []
    for col, pair in filters.items():
        low, high = pair
        vals = mads_df[col]
        if low is not None and high is not None:
            masks.append((vals > high) | (vals < low))
        elif low is not None:
            masks.append(vals < low)
        elif high is not None:
            masks.append(vals > high)
    if reduction == "all":
        is_outlier = reduce(lambda x, y: x & y, masks)
    else:
        is_outlier = reduce(lambda x, y: x | y, masks)
    filtered = adata[~is_outlier, :]
    failed = adata[is_outlier]
    return filtered, failed
