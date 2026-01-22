#!/usr/bin/env ipython

from collections.abc import Sequence
from pathlib import Path
from typing import Literal

import anndata as ad
import matplotlib.pyplot as plt
import mudata as md
import numpy as np
import pandas as pd
import plotnine as gg
import polars as pl
import polars.selectors as cs
import scirpy as ir
import seaborn as sns
from matplotlib.figure import Figure
from plotnine.ggplot import ggplot

from alignment import align_vdj, get_stitchr_seqs, plot_vdj
from analyses import SCOL, get_airr, maybe_filter_by_rank

md.set_options(pull_on_update=False)


try:
    from snakemake.script import snakemake as smk
except ImportError:
    smk = type("snakemake", (), {"rule": None, "config": {}})

SCOL: str = smk.config.get("sample_col", "Sample_Name")


def categorical_to_str(
    df: pd.DataFrame, subset: Sequence | None = None
) -> pd.DataFrame:
    cols = subset or df.columns
    new_df = df.copy()
    to_cast = {c: str for c in cols if hasattr(new_df[c], "cat")}
    new_df = new_df.astype(to_cast)
    for c in to_cast:
        if any(new_df[c].str.contains("^[Nn]a[Nn]?$")):
            new_df[c] = new_df[c].replace("^[Nn]a[Nn]?$", np.nan, regex=True)
    return new_df


def from_pandas_categorical(
    df: pd.DataFrame, mapping: dict | None = None, make_categorical: bool = True
) -> pl.DataFrame:
    "Convert a pandas df to to polars df while handling categorical types"
    mapping = mapping or {}
    type_series = df.dtypes
    categoricals = type_series[
        [isinstance(x, pd.CategoricalDtype) for x in type_series]
    ].index
    type_mapping = {k: mapping.get(k, "str") for k in categoricals}
    converted = pl.from_pandas(df.astype(type_mapping))
    if not make_categorical:
        return converted
    return converted.cast(
        {k: pl.Categorical for k in type_mapping.keys() if k not in mapping}
    )


def plot_clone_ranking(
    mdata: md.MuData,
    k: int = 5,
    expanded_in: str | None = None,
    target_col: str = "clone_id",
    airr_mod: str = "airr",
) -> ggplot:
    """Plot the proportion of cells made up by the top k clones in the repertoire

    Parameters
    ----------
    expanded_in : str | None
        Calculate clone ranks with respect to this column e.g. specify samples to
        show top clones for each sample
    k : int
        Top k clones to include in the plot

    """
    to_select = [target_col] if not expanded_in else [target_col, expanded_in]
    df = from_pandas_categorical(mdata[airr_mod].obs.loc[:, to_select])
    top_clones = pl.concat(
        [
            d[1][target_col]
            .value_counts()
            .with_columns(
                pl.col("count").rank("dense", descending=True).alias("rank"),
                (pl.col("count") / pl.col("count").sum()).alias("prop"),
            )
            .with_columns(
                pl.when(pl.col("rank") > k)
                .then(pl.lit("other"))
                .otherwise("rank")
                .alias("rank")
            )
            .group_by("rank")
            .agg(
                pl.col("prop").sum(),
                pl.len().alias("rank_count"),
                pl.col("count").sum(),
            )
            .sort("rank")
            .with_columns(pl.lit(d[0][0]).alias(expanded_in))
            for d in df.group_by(expanded_in)
        ]
    ).with_columns(
        pl.when(pl.col("rank_count") == 1)
        .then(pl.lit(None))
        .otherwise("rank_count")
        .cast(pl.String)
        .alias("rank_count")
    )
    kws = {"fill": "rank", "x": expanded_in, "y": "prop", "label": "rank_count"}
    plot = (
        gg.ggplot(top_clones, gg.aes(**kws))
        + gg.geom_bar(position="stack", stat="identity")
        + gg.geom_text(position=gg.position_fill(vjust=0.5))
        + gg.labs(fill="Clone Rank")
        + gg.ylab("Proportion")
        + gg.ggtitle(
            title=f"Proportion of top {k} clones",
            subtitle="Numbers within bars are the count of unique clones in each rank",
        )
    )
    return plot


def plot_alpha_diversity(
    data: ad.AnnData | md.MuData, metrics, groupby: str, airr_mod: str = "airr"
) -> ggplot:
    """
    Plot multiple diversity metrics. `data` is assumed to have had metrics pre-computed
    Works best if metrics are normalized to the same scale
    """
    if isinstance(data, md.MuData):
        data = data[airr_mod]
    to_rename = {f"{x}_clone_id": x for x in metrics}
    df = (
        data.obs.loc[:, [groupby] + list(to_rename.keys())]
        .rename(to_rename, axis=1)
        .drop_duplicates()
        .melt(id_vars=groupby)
    )
    plot = (
        gg.ggplot(df, gg.aes(x=groupby, y="value", fill="variable"))
        + gg.geom_bar(position=gg.position_dodge(), stat=gg.stat_identity())
        + gg.labs(fill="Metric")
        + gg.ylab("Value")
    )
    return plot


def plot_group_abundance(
    data: ad.AnnData | md.MuData,
    x: str,
    ylab: str = "Number of Cells",
    xlab: str | None = None,
    fill: str = "clonal_expansion",
    normalize: bool = False,
    max_cols: int = 50,
    at_least: int = 20,
    ignore_nan: bool = True,
):
    "Custom plot because the scirpy version won't show legends for some reason"
    data.obs = categorical_to_str(data.obs, [fill, x])
    xlab = xlab if xlab is not None else x
    n_cols: pd.Series = data.obs[x].value_counts()
    if len(n_cols) > max_cols:
        filtered: pd.DataFrame = data.obs.loc[
            data.obs[x].isin(n_cols[:max_cols].index), :
        ]
    else:
        filtered = data.obs
    if ignore_nan:
        mask = filtered[fill].isna() | filtered[x].isna()
        filtered = filtered.loc[~mask, :]
    if ":" in x:
        filtered = filtered.rename({x: x.replace(":", "_")}, axis=1)
        x = x.replace(":", "_")
    if at_least > 0:
        counts = pd.crosstab(filtered[x], filtered[fill]).sum(axis=1)
        passed_vals = set(counts[counts >= at_least].index)
        filtered = filtered.loc[filtered[x].isin(passed_vals), :]
    position = "fill" if normalize else "stack"
    plot = (
        gg.ggplot(
            filtered,
            gg.aes(x=f"reorder({x}, {x}, len)", fill=fill),
        )
        + gg.geom_bar(position=position)
        + gg.ylab(ylab)
        + gg.xlab(xlab)
    )
    return plot


def plot_obs(
    adata: ad.AnnData,
    mode: Literal["tSNE", "UMAP"],
    hues: Sequence,
    from_bd_pipeline: bool = True,
    **kwargs,
) -> Figure:
    fig, ax = plt.subplots(1, len(hues))
    if mode == "tSNE" and from_bd_pipeline:
        x = "tSNE_1"
        y = "tSNE_2"
    elif mode == "UMAP" and from_bd_pipeline:
        x = "UMAP_1"
        y = "UMAP_2"
    else:
        raise NotImplementedError()
    for i, hue in enumerate(hues):
        sns.scatterplot(data=adata.obs, x=x, y=y, hue=hue, ax=ax[i], **kwargs)
        ax[i].set_ylabel(x.replace("_", ""))
        ax[i].set_xlabel(y.replace("_", ""))
        ax[i].get_yaxis().set_ticks([])
        ax[i].get_xaxis().set_ticks([])
        if i > 0:
            ax[i].set_ylabel(None)
    return fig


def top_clone_calls(
    adata: ad.AnnData, k: int = 5, key: str = "clone_id"
) -> pd.DataFrame:
    """Return a dataframe containing the V, D, J gene calls for
        the most abundant clonotypes or clone clusters

    Parameters
    ----------
    k : int
        Return results for the top k clonotypes
    adata : ad.AnnData
        AnnData object with IR
    """
    cols = ["sequence_id", "c_call", "d_call", "j_call", "v_call"]
    with ir.get.airr_context(adata, cols) as m:
        top = adata.obs[key].value_counts()[:k].reset_index()
        top["rank"] = top["count"].rank(ascending=False, method="dense")
        call_cols = [col for col in m.obs.columns if "_call" in col]
        filtered = m.obs.loc[m.obs[key].isin(top[key]), :].loc[:, [key] + call_cols]
        df = top.merge(filtered, on=key).drop_duplicates().astype({"rank": int})
    return df


# * Generate plots
if smk.rule == "make_reports":
    mdata: md.MuData = md.read_h5mu(smk.input[0])
    if to_ignore := smk.config.get("ignore_samples"):
        mdata = mdata[~mdata.obs[SCOL].isin(to_ignore), :]
    ct_col = smk.config["cell_type_col"]
    airrs = {
        s: mdata[mdata.obs[SCOL] == s, :]
        for s in mdata.obs[SCOL].unique()
        if s not in to_ignore
    }
    # TODO: Plot on per-run basis
    # tsne_fig = plot_obs(adata, "tSNE", hues=to_plot)
    # umap_fig = plot_obs(adata, "UMAP", hues=to_plot)
    clone_calls, cluster_calls = [], []
    for p in ["expansion", "ct_abundance"]:
        if not (check_path := Path(smk.output[p])).exists():
            check_path.mkdir(parents=True)

    for sample, cur in airrs.items():
        clone_expansion_plot = (
            plot_group_abundance(
                cur, x=ct_col, fill="airr:clonal_expansion", max_cols=30, at_least=10
            )
            + gg.guides(fill=gg.guide_legend(title="Clone Size", reverse=True))
            + gg.theme(axis_text_x=gg.element_text(rotation=45, hjust=1))
            + gg.ggtitle(sample)
        )
        clone_expansion_plot.save(
            f"{smk.output['expansion']}/{sample}.png",
            width=10,
            height=10,
            verbose=False,
        )

        ct_plot = (
            plot_group_abundance(
                cur, x="airr:clone_id", fill=ct_col, max_cols=20, at_least=2
            )
            + gg.xlab("Clone ID")
            + gg.ggtitle(sample)
        )
        ct_plot.save(
            f"{smk.output['ct_abundance']}/{sample}.png",
            width=13,
            height=10,
            verbose=False,
        )
        for key, lst in zip(["clone_id", "cc_id"], [clone_calls, cluster_calls]):
            tmp = top_clone_calls(cur["airr"], k=5, key=key).assign(**{SCOL: sample})
            lst.append(tmp.iloc[:, [-1] + list(range(len(tmp.columns) - 2))])

    c_rank_plot = plot_clone_ranking(mdata, k=4, expanded_in=SCOL)
    c_rank_plot.save(smk.output["clone_ranks"], width=10, height=8, verbose=False)

    cc_rank_plot = plot_clone_ranking(
        mdata, k=4, expanded_in=SCOL, target_col="cc_id"
    ) + gg.labs(fill="Clone Cluster Rank")
    cc_rank_plot.save(
        smk.output["clone_cluster_ranks"], width=10, height=8, verbose=False
    )

    dplot = plot_alpha_diversity(
        mdata["airr"], smk.config["alpha_metrics"], groupby=SCOL
    ) + gg.xlab("Sample")
    dplot.save(smk.output["alpha_diversity"], verbose=False)
    public_private_plot = plot_group_abundance(
        mdata["airr"], x="clone_id", fill=SCOL, max_cols=30
    ) + gg.xlab("Clone ID")
    public_private_cc_plot = plot_group_abundance(
        mdata["airr"], x="cc_id", fill=SCOL, max_cols=30
    ) + gg.xlab("Cluster ID")

    public_private_plot.save(
        smk.output["public_private"], width=15, height=8, verbose=False
    )
    public_private_cc_plot.save(
        smk.output["public_private_clusters"], width=15, height=8, verbose=False
    )
    pd.concat(clone_calls).to_csv(smk.output["top_clones"], index=False)
    pd.concat(cluster_calls).to_csv(smk.output["top_clusters"], index=False)
if smk.rule == "plot_sequence":
    outdir = Path(smk.params["outdir"])
    plotdir = Path(smk.params["plotdir"])
    airr = get_airr(filter_samples=True)
    viz_config = smk.config["vdj_plot"]
    for_title = [c[1] for c in viz_config["title_spec"]]
    if "clone_id" not in for_title:
        for_title.append("clone_id")
    _, airr = maybe_filter_by_rank(airr, viz_config)
    cols = ["sequence", "cdr3", "cdr1", "cdr2", "fwr1", "fwr2", "fwr3", "fwr4"] + [
        s
        for seq in [
            [f"{k}_sequence_start", f"{k}_sequence_end", f"{k}_call"]
            for k in ("v", "d", "j")
        ]
        for s in seq
    ]
    cols = ["sequence_id"] + cols + ["c_call"]
    allele_df = pl.from_pandas(
        ir.get.airr(airr, ["v_call", "c_call"]).reset_index(names="index")
    ).with_columns(cs.ends_with("c_call") + "*01")
    # WARNING: BD pipeline doesn't provide allele specificity for C, so we will just
    # take the first allele
    cohort_v_alleles = list(allele_df["VJ_1_v_call"]) + list(allele_df["VDJ_1_v_call"])
    cohort_c_genes = list(allele_df["VJ_1_c_call"]) + list(allele_df["VDJ_1_c_call"])
    v_leaders: pl.DataFrame = get_stitchr_seqs(
        set(cohort_v_alleles), seqtype="~LEADER"
    ).rename({"sequence": "v_leader"})
    c_genes: pl.DataFrame = get_stitchr_seqs(
        set(cohort_c_genes), seqtype="~CONSTANT"
    ).rename({"sequence": "c_gene"})
    c_gene_dict = dict(zip(c_genes["allele"], c_genes["c_gene"]))
    c_genes = (
        allele_df.select(pl.col("index"), cs.ends_with("c_call"))
        .with_columns(
            cs.ends_with("c_call").map_elements(
                lambda x: c_gene_dict.get(x, None), return_dtype=pl.String
            )
        )
        .rename(
            lambda x: x if "c_call" not in x else x.replace("_c_call", "_c_sequence")
        )
    )
    for sample in airr.obs[SCOL].unique():
        mask = airr.obs[SCOL] == sample
        cur = airr[mask, :]
        seqs = pl.from_pandas(ir.get.airr(cur, cols).reset_index(names="index"))
        title_obs = pl.from_pandas(cur.obs.loc[:, for_title])
        seqs = pl.concat([seqs, title_obs], how="horizontal")
        for chain in ("VJ_1", "VDJ_1"):
            id_col = f"{chain}_sequence_id"
            chain_seqs: pl.DataFrame = (
                seqs.select(cs.starts_with(chain), pl.col(for_title + ["index"]))
                .unique([f"{chain}_j_call", f"{chain}_v_call", f"{chain}_cdr3"])
                .with_columns(
                    pl.struct(["clone_id", id_col])
                    .map_elements(
                        lambda x: f"cid{x['clone_id']}-{x[id_col]}",
                        return_dtype=pl.String,
                    )
                    .alias(id_col)
                )
            )
            chain_seqs = chain_seqs.join(
                v_leaders, left_on=f"{chain}_v_call", right_on="allele", how="left"
            )
            cur_c = c_genes.select(cs.starts_with(chain), pl.col("index"))
            chain_seqs = (
                chain_seqs.join(cur_c, on="index", how="left")
                .filter(pl.col(id_col).is_not_null())
                .with_columns(
                    pl.col(f"{chain}_c_sequence").str.head(
                        pl.col(f"{chain}_sequence").str.len_chars()
                        - pl.col(f"{chain}_j_sequence_end")
                    )
                )
            )
            cur_outdir = outdir / sample / chain
            cur_plotdir = plotdir / sample / chain
            cur_plotdir.mkdir(exist_ok=True, parents=True)
            successes = align_vdj(
                chain_seqs,
                chain,
                outdir=cur_outdir,
                id_col=id_col,
                additional_seqs={
                    "v_leader": "v_leader",
                    "c_gene": f"{chain}_c_sequence",
                },
            )
            for afile in cur_outdir.glob("*.fasta"):
                outfile = cur_plotdir / f"{afile.stem}.png"
                fig = plot_vdj(
                    afile.stem,
                    chain,
                    chain_seqs,
                    file=afile,
                    id_col=id_col,
                    wrap_length=viz_config.get("wrap_length", 200),
                    title_spec=[("", id_col)]
                    + viz_config["title_spec"]
                    + [("C call", f"{chain}_c_call")],
                    color_scheme="Clustal",
                )
                fig.set_figwidth(viz_config.get("width", 30))
                fig.set_dpi(viz_config.get("dpi", 100))
                fig.savefig(outfile)
