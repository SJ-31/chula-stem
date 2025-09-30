#!/usr/bin/env ipython

from collections.abc import Sequence
from typing import Literal

import anndata as ad
import matplotlib.pyplot as plt
import mudata as md
import pandas as pd
import plotnine as gg
import seaborn as sns
from matplotlib.figure import Figure
from plotnine.ggplot import ggplot
from pynndescent.sparse import diversify

try:
    from snakemake.script import snakemake as smk
except ImportError:
    smk = type("snakemake", (), {"rule": None})


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
    ignore_nan: bool = True,
):
    "Custom plot because the scirpy version won't show legends for some reason"
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
    plot = (
        gg.ggplot(
            filtered,
            gg.aes(x=f"reorder({x}, {x}, len)", fill=fill),
        )
        + gg.geom_bar()
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


# * Generate plots
if smk.rule == "make_plots":
    mdata: md.MuData = md.read_h5mu(smk.input[0])
    ct_col = smk.config["ct_col"]
    airrs = {
        s: mdata[mdata.obs["Sample_Name"] == s, :]
        for s in mdata.obs["Sample_Name"].unique()
        if s not in {"Undetermined", "Multiplet"}
    }
    # TODO: Plot on per-run basis
    # tsne_fig = plot_obs(adata, "tSNE", hues=to_plot)
    # umap_fig = plot_obs(adata, "UMAP", hues=to_plot)

    for sample, cur in airrs.items():
        clone_expansion_plot = (
            plot_group_abundance(cur, x=ct_col, fill="airr:clonal_expansion")
            + gg.theme_classic()
            + gg.guides(fill=gg.guide_legend(title="Clone Size", reverse=True))
        )
        clone_expansion_plot.save()  # TODO:

        ct_plot = plot_group_abundance(
            cur, x="airr:clone_id", fill=ct_col, max_cols=20
        ) + gg.xlab("Clone ID")
        ct_plot.save()  # TODO:

    dplot = plot_alpha_diversity(
        mdata, smk.config["alpha_metrics"], groupby="Sample_Name"
    )
    dplot.save()
    public_private_plot = plot_group_abundance(
        mdata, x="airr:clone_id", fill="Sample_Name", max_cols=20
    )
    public_private_plot.save()  # TODO:
