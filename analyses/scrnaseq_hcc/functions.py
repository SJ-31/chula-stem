#!/usr/bin/env python

import re
from functools import reduce
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Callable, Literal

import anndata as ad
import decoupler as dc
import marimo as mo
import numpy as np
import pandas as pd
import plotnine as gg
import polars as pl
import pymupdf
import scanpy as sc
import snakemake.io
import yaml
from chula_stem.r_utils import edgeR_ovr
from chula_stem.sc_rnaseq import annotate_adata_vars, annotate_marker, distance_by_mads
from chula_stem.utils import read_existing
from pymupdf import Document

# * Data prep/retrieval


def prepare_data(file, feature_file, env):
    droot = Path(env["data_root"])
    tmp = []
    manifest = pl.read_csv(env["manifest"])
    target = "filtered_feature_bc_matrix.h5"
    for row in manifest.iter_rows(named=True):
        patient = row["sample_name"]
        stype = row["type"]
        suffix = f"{row['type']}-{row['treatment']}"
        if not (file_key := env.get("manifest_file_key")):
            data_path: Path = (
                droot / patient / "processed" / f"cellranger_{suffix}/{target}"
            )
            if not data_path.exists():
                print(
                    "The files for sample {} ({}) don't exist".format(patient, suffix)
                )
                continue
            current: ad.AnnData = sc.read_10x_h5(data_path)
        else:
            current = ad.read_h5ad(row[file_key])
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
                print(f"WARNING: Path to extra file {path} does not exist")
                continue
            current = ad.read_h5ad(path)
            current.obs.loc[:, "patient"] = patient
            current.obs.loc[:, "sample"] = f"{patient}_{stype}-{treatment}"
            current.obs.loc[:, "type"] = stype
            current.obs.loc[:, "treatment"] = treatment
            current.obs.loc[:, "flowcell"] = "unknown"
            tmp.append(current)
    adata = ad.concat(
        tmp, merge="first", index_unique="-", join="outer", uns_merge="same"
    )
    adata = adata[:, ~adata.var["gene_ids"].isna()]
    annotate_adata_vars(adata, "gene_ids", savepath=Path(env["gene_reference"]))
    marker_genes = env.get("obs_markers_annotate")
    if marker_genes:
        annotate_marker(adata, marker_genes=marker_genes, gene_col="hgnc_symbol")
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
    hvgs = permissive_feature_selection(adata, env)
    sc.pp.pca(
        adata,
        n_comps=100,
        layer="lshift_normalized",
        mask_var=hvgs["mask"],
    )
    with open(feature_file, "w") as f:
        f.write("\n".join(hvgs["index"]))
    adata.write_h5ad(file)
    adata = ad.read_h5ad(file, backed=True)
    return adata


def fs_dispatch(method: str, adata: ad.AnnData) -> Callable:
    if method == "seurat":
        return lambda **kws: sc.pp.highly_variable_genes(
            adata=adata, flavor="seurat", **kws
        )
    raise NotImplementedError()


def permissive_feature_selection(adata: ad.AnnData, env: dict) -> dict:
    cfg = env["permissive_fs"]
    kws = cfg.get("kws") or {}
    result = {}
    method = cfg["method"]
    selection_fn = fs_dispatch(method, adata)
    features = selection_fn(**kws)
    if method in {"seurat", "cellranger"}:
        result["mask"] = features["highly_variable"].values
        result["index"] = features.query("highly_variable").index.tolist()
    else:
        raise NotImplementedError()
    return result


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
        tmp = pd.read_csv(old_clusters).rename(columns={"cluster": "tmp_cluster"})
        cells2old_clusters = pd.Series(tmp["tmp_cluster"])
        cells2old_clusters.index = tmp["cell_id"]
        old = cells2old_clusters[adata.obs_names.tolist()]
    else:
        old = adata.obs[old_clusters]
    communities = pd.read_csv(cfs_community_file)
    lookup = {(s, g): c for _, (s, g, c) in communities.iterrows()}
    adata.obs = adata.obs.merge(old, left_index=True, right_index=True)
    reassigned = adata.obs["sample"].combine(
        adata.obs["tmp_cluster"], lambda x, y: lookup.get((x, y))
    )
    adata.obs.loc[:, column] = reassigned.astype(str)


def get_gs_and_cell_markers(
    env: dict,
    min_count: tuple[int | None, int | None] | None = None,
    as_df: bool = False,
    df_names=("set", "symbol"),
) -> tuple[dict | pd.DataFrame, dict | pd.DataFrame]:
    cell_markers = env["markers"] or {}
    gene_sets = env["gene_sets"] or {}
    min_count = min_count or (None, None)
    for update_to, file_list, minimum in zip(
        (gene_sets, cell_markers), ("gene_set_files", "marker_files"), min_count
    ):
        for file in env[file_list]:
            with open(file, "r") as f:
                tmp = yaml.safe_load(f)
                if minimum is not None:
                    tmp = {
                        k: v
                        for k, v in tmp.items()
                        if (isinstance(v, list) and len(v) > minimum)
                        or (minimum <= 1 and isinstance(v, str))
                    }
                update_to.update(tmp)
    if as_df:
        result = [
            pd.DataFrame({df_names[0]: dct.keys(), df_names[1]: dct.values()})
            .explode(df_names[1])
            .drop_duplicates()
            for dct in [gene_sets, cell_markers]
        ]
        gene_sets, cell_markers = result
    return gene_sets, cell_markers


def with_normalized_layers(adata: ad.AnnData, env: dict) -> None:
    "Add additional normalization layers to `adata`, retaining the raw data in X"
    cfg = env["normalization"]
    # BUG: this wasn't working
    lshift = sc.pp.normalize_total(adata, inplace=False, **cfg["normalize_total"])
    adata.layers["lshift_normalized"] = lshift["X"]
    adata.obs["lshift_size_factors"] = lshift["norm_factor"]
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


# ** Specifying output


def get_scvi_model_path(fs_name: str, prefix: str, env: dict) -> dict:
    return {
        "dir_path": f"{env['outdir']}/{fs_name}/.scvi",
        "prefix": f"{prefix}_",
        "overwrite": True,
    }


def provide_output_from_fs(fs_name: str, env: dict) -> dict:
    """
    Specify all output paths downstream of feature (highly variable gene) selection

    The keys are always prefixed with FEATURE_SELECTION_METHOD_NAME-
    """
    root: Path = Path(env["outdir"]) / fs_name
    prefix = f"{fs_name}-"
    outs = {}
    integration_methods = list(env["integration"].keys())
    integration_methods.append("unintegrated")
    # Integration methods
    outs[f"{prefix}_dotplots_PLOT"] = f"{root}/dotplots.pdf"
    outs[f"{prefix}_cellassign_PLOT"] = f"{root}/cellassign_metrics.pdf"

    # DR results
    # Key format is FEATURE_SELECTION-INTEGRATION_dr_DR_NAME
    for imethod in integration_methods:
        if imethod != "unintegrated":
            outs[f"{prefix}{imethod}"] = f"{root}/{imethod}_integrated.h5ad"
            outs[f"{prefix}{imethod}_clustering"] = f"{root}/{imethod}_clustering.h5ad"
            if not env["test"]:
                outs[f"{prefix}{imethod}_clustering_PLOT"] = (
                    f"{root}/{imethod}_clustering.pdf"
                )

        for method, values in env["DR"]["methods"].items():
            v = values["vary"][1]
            outs[f"{prefix}{imethod}_dr_{method}"] = snakemake.io.expand(
                f"{root}/DR_{imethod}/{method}/{method}_{{v}}.npy", v=v
            )
    # CFS (hvgs is done independently for each sample here)
    for file in ("similarity", "communities", "ind_clustering"):
        outs[f"{prefix}{file}"] = f"{root}/cfs/{file}.csv"

    # DE analysis
    # TODO: specify this

    return outs


def provide_annotation_output(env) -> dict:
    outdir = env["outdir"]
    root = f"{outdir}/annotations"
    result = {
        "cellassign_predictions": f"{outdir}/cellassign_predictions.csv",
        "cellassign_run_metrics": f"{outdir}/cellassign_metrics.csv",
    }
    if not env.get("chosen_clusters"):
        return result
    for csv in (
        "gene_set_activity",
        "marker_gene_activity",
        "clusters-edgeR_de",
        "clusters-scVI_de",
        "samples_de.csv",
    ):
        if csv.startswith("clusters") and not env.get("chosen_clusters"):
            continue
        elif csv.startswith("samples_de") and not env.get("do_de_samples"):
            continue
        result[csv] = f"{root}/{csv}.csv"
    return result


# * Plotting


def plot_clusters_in_samples(
    adata, cluster_col: str, sample_col: str = "sample", ncol: int = 3
) -> gg.ggplot:
    """Visualize the size of each cluster within each sample"""
    adata = adata[~adata.obs[cluster_col].isna(), :]
    cells_per_sample = adata.obs[sample_col].value_counts()
    bulked = dc.pp.pseudobulk(
        adata, sample_col=sample_col, groups_col=cluster_col, layer=None, mode="sum"
    )
    bulked.obs.loc[:, "psbulk_cells_prop"] = (
        bulked.obs["psbulk_cells"].values
        / cells_per_sample[bulked.obs[sample_col]].values
    )
    plot = (
        gg.ggplot(
            bulked.obs,
            gg.aes(
                x=cluster_col,
                y="psbulk_cells_prop",
                fill="psbulk_counts",
            ),
        )
        + gg.geom_bar(stat="identity")
        + gg.facet_wrap(sample_col, drop=False, dir="v", ncol=ncol)
        + gg.guides(
            fill=gg.guide_legend(title="Total count"),
        )
        + gg.ylab("Proportion in sample")
        + gg.theme_bw()
        + gg.theme(
            panel_grid_major=gg.element_blank(), panel_grid_minor=gg.element_blank()
        )
    )
    return plot


def make_cluster_dotplots(
    adata: ad.AnnData,
    filename: str | Path,
    markers: dict | list,
    env: dict,
    cluster_results: str | Path | None = None,
    group_rotation: int = 0,
    additional_groups: list | None = None,
    with_samples: bool = True,
) -> None:
    """Save a dotplot for each clustering sweep for a specified combination of
    (feature selection method, integration_method) to a single file

    Parameters
    ----------
    cluster_results : str | Path
        Path to a h5ad file containing the clustering results in obs. Must have the same
        index as adata
    markers : dict | list
        vars to plot
    additional_groups : list | None
        Additional groupings to show on dotplots e.g. cellassign labels
    group_rotation : int
        Rotation for group labels when on the x axis
    """
    if cluster_results:
        loaded: ad.AnnData = ad.read_h5ad(cluster_results, backed=True)
        adata.obs = adata.obs.merge(
            loaded.obs, left_index=True, right_index=True, how="left"
        )
        group_cols = loaded.obs.columns
    else:
        group_cols = []
    if groups_from_env := env.get("additional_groups"):  # Can put cellassign in here
        additional_groups = additional_groups + groups_from_env
    if additional_groups:
        group_cols = additional_groups + group_cols
    if len(group_cols) == 0:
        raise ValueError("No grouping columns could be specified")
    cfg: dict = env.get("dotplot") or {}
    kws = cfg.get("kws") or {}
    totals_kws = cfg.get("totals_kws") or {}
    if isinstance(markers, dict) and (from_cfg := cfg.get("markers")):
        print(f"INFO: Adding additional markers from dotplot configuration: {from_cfg}")
        markers.update(from_cfg)
    doc: Document = pymupdf.open()
    transpose = kws.get("swap_axes", False)
    with TemporaryDirectory() as tmp:
        for i, col in enumerate(group_cols):
            save_to_1 = f"{tmp}/{i}_1.pdf"
            save_to_2 = f"{tmp}/{i}_2.pdf"
            call = {
                "adata": adata[adata.obs[col].notnull(), :],
                "title": f"Groups: {col}",
                "groupby": col,
                "var_names": markers,
                "show": False,
                "return_fig": True,
                **kws,
            }
            plot = sc.pl.dotplot(**call).add_totals(**totals_kws)
            plot.make_figure()
            if transpose:
                axes: dict = plot.get_axes()
                main = axes["mainplot_ax"]
                for label in main.get_xticklabels():
                    label.set_rotation(group_rotation)
                # BUG: the rotation works, but you need to shift the text a bit and
                # and there's an issue that scanpy will automatically truncate the labels
                # if gene_groups := axes.get("gene_group_ax"):
                #     for child in gene_groups.get_children():
                #         if isinstance(child, Text):
                #             if child.get_text():
                #                 child.set_rotation(0)
            plot.fig.savefig(save_to_1, bbox_inches="tight")
            doc.insert_file(save_to_1)
            if col != "sample" and with_samples:
                cluster_counts = plot_clusters_in_samples(
                    adata, col, ncol=2
                ) + gg.theme(
                    figure_size=(15, 10), axis_text_x=gg.element_text(rotation=90)
                )
                cluster_counts.save(save_to_2)
                doc.insert_file(save_to_2)
    doc.save(filename)


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


def add_saved_dr(
    adata: ad.AnnData, env: dict, fs_name: str, integration: str = "unintegrated"
) -> None:
    outs = provide_output_from_fs(fs_name, env)
    for method in env["DR"]["methods"].keys():
        key = f"{fs_name}-{integration}_dr_{method}"
        files = outs[key]
        if files:
            directory = Path(files[0]).parent
            if directory.exists():
                for efile in directory.iterdir():
                    if efile.suffix == ".npy":
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


# * Sample-level DE


# * Cluster analysis


def add_clusterings(adata, clusters_to_add: dict[str, list], env: dict) -> list:
    outdir = Path(env["outdir"])
    previously_read: dict = {}
    added = []
    for name, (fs_name, integration, cluster_name) in clusters_to_add.items():
        if integration is None and cluster_name == "cfs":
            add_cfs_community(
                adata,
                old_clusters=outdir / fs_name / "cfs/ind_clustering.csv",
                cfs_community_file=outdir / fs_name / "cfs/communities.csv",
                column=name,
            )
            added.append(name)
            continue
        cluster_results: pd.DataFrame | None = previously_read.get(
            (fs_name, integration)
        )
        if cluster_results is None:
            path = outdir / fs_name / f"{integration}_clustering.h5ad"
            obs = ad.read_h5ad(path).obs
            cluster_results = adata.obs.merge(
                obs, left_index=True, right_index=True, how="left"
            )
            previously_read[(fs_name, integration)] = cluster_results
        adata.obs.loc[:, name] = cluster_results[cluster_name]
        added.append(name)
    return added


# ** Enrichment


def enrich_clusters(
    adata,
    cfg: dict,
    env: dict,
    gs_out: str,
    marker_out: str,
):
    """Determine which clusters are enriched in the gene sets
        specified in `gene_set_files` and `gene_sets`

    Notes
    -----
    Routine steps
    1. Compute per-cell enrichment scores for the specified gene sets
    2. Use dc to get a consensus score
    3. For each cluster, rank the activity of each gene set using t-test against
    the the activity in all other clusters combined

    Test statistic by default is T statistic, (mean1-mean2)/se,
    Which you should be able to safely interpret as the magnitude of the score
    They use two-sided test, so convert to absolute
    Positive means gene set under `name` has high expression and vice-versa

    """
    gs_df, marker_df = get_gs_and_cell_markers(env, None, True, ("source", "target"))
    if exclude := env["cluster_cells"].get("exclude"):
        adata = adata[~adata.obs["sample"].isin(exclude), :]
    cluster_names = add_clusterings(
        adata, clusters_to_add=env["chosen_clusters"], env=env
    )
    kws = cfg.get("decouple_kws") or {}
    cutoff = cfg.get("cutoff", 0.0001)
    for out, df in zip((gs_out, marker_out), (gs_df, marker_df)):
        ranked_dfs = []
        dc.mt.decouple(adata, df, **kws)
        methods = kws["methods"]
        if len(methods) > 1:
            dc.mt.consensus(adata)
            key = "score_consensus"
        else:
            key = f"score_{methods[0]}"
        scores: ad.AnnData = dc.pp.get_obsm(adata, key=key)
        for clst in cluster_names:
            ranked = dc.tl.rankby_group(scores, groupby=clst)
            ranked = ranked.loc[ranked["pval"] <= cutoff, :].assign(clustering=clst)
            ranked_dfs.append(ranked)
        pd.concat(ranked_dfs).to_csv(out, index=False)


# ** DE


def do_de_clusters(method: str, adata: ad.AnnData, cfg: dict, env: dict, **kws) -> dict:
    if exclude := env["cluster_cells"].get("exclude"):
        adata = adata[~adata.obs["sample"].isin(exclude), :]
    cluster_names = add_clusterings(
        adata, clusters_to_add=env["chosen_clusters"], env=env
    )
    extra_contrasts = cfg.get("extra_contrasts")
    cluster_kws = cfg.get("kws") or {}
    if method == "edgeR":
        return de_clusters_edgeR(adata, cluster_names, extra_contrasts, cluster_kws)
    elif method == "scVI":
        return de_clusters_scVI(
            adata, cluster_names, extra_contrasts, cluster_kws, **kws
        )
    raise ValueError(f"unsupported method {method}")


def de_clusters_scVI(
    adata: ad.AnnData,
    cluster_names,
    extra_contrasts: list | None,
    kws: dict,
    model_file,
    feature_file,
) -> dict:
    import scvi

    with open(feature_file, "r") as f:
        features = f.read().splitlines()
    adata = adata[:, features]
    model: scvi.model.SCVI = scvi.model.SCVI.load(model_file, adata)
    top_de = []
    for clst in cluster_names:
        result = model.differential_expression(groupby="clst", **kws)
        top_de.append(result.assign(contrast=clst))
    if extra_contrasts:
        for contrast in extra_contrasts:
            result = model.differential_expression(**contrast, **kws)
            top_de.append(result.assign(contrast=contrast))
    return {"top_de": pd.concat(top_de)}


def de_clusters_edgeR(adata: ad.AnnData, cluster_names, extra_contrasts, kws) -> dict:
    """Perform DE analysis between clusters (OVR) using edgeR"""
    de_counts = []
    all_top = []
    for clst in cluster_names:
        bulked = dc.pp.pseudobulk(
            adata[~adata.obs[clst].isna()], "sample", groups_col=clst
        )
        sc.pp.filter_cells(bulked, min_counts=15)
        sc.pp.filter_genes(bulked, min_counts=50)
        n_de, top_de = edgeR_ovr(
            bulked, clst, extra_contrasts=extra_contrasts.get(clst), **kws
        )
        de_counts.append(
            n_de.reset_index(names="direction")
            .melt("direction")
            .rename({"variable": "cluster", "value": "count"})
            .assign(clustering=clst)
        )
        all_top.append(top_de.assign(clustering=clst))
    return {"de_counts": pd.concat(de_counts), "top_de": pd.concat(all_top)}


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
