#!/usr/bin/env ipython

from collections.abc import Sequence
from typing import Literal

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import plotnine as gg
import seaborn as sns
from matplotlib.figure import Figure

try:
    from snakemake.script import snakemake as smk
except ImportError:
    smk = {}


def plot_group_abundance(
    adata: ad.AnnData,
    x: str,
    ylab: str = "Number of Cells",
    xlab: str | None = None,
    fill: str = "clonal_expansion",
    normalize: bool = False,
    max_cols: int = 50,
):
    "Custom plot because the scirpy version won't show legends for some reason"
    xlab = xlab if xlab is not None else x
    n_cols: pd.Series = adata.obs[x].value_counts()
    if len(n_cols) > max_cols:
        filtered: pd.DataFrame = adata.obs.loc[
            adata.obs[x].isin(n_cols[:max_cols].index), :
        ]
    else:
        filtered = adata.obs
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
    adata = ad.read_h5ad(smk.input[0])
    ct_col = smk.config["ct_col"]
    airrs = {
        s: adata[adata.obs["Sample_Name"] == s, :]
        for s in adata.obs["Sample_Name"].unique()
        if s not in {"Undetermined", "Multiplet"}
    }
    # TODO: Plot on per-run basis
    # tsne_fig = plot_obs(adata, "tSNE", hues=to_plot)
    # umap_fig = plot_obs(adata, "UMAP", hues=to_plot)

    for sample, cur in airrs.items():
        clone_expansion_plot = (
            plot_group_abundance(cur, x=ct_col, fill="clonal_expansion")
            + gg.theme_classic()
            + gg.guides(fill=gg.guide_legend(title="Clone Size", reverse=True))
        )
        clone_expansion_plot.save()  # TODO:

        ct_plot = plot_group_abundance(
            cur, x="clone_id", fill=ct_col, max_cols=20
        ) + gg.xlab("Clone ID")
        ct_plot.save()  # TODO:

    public_private_plot = plot_group_abundance(
        adata, x="clone_id", fill="Sample_Name", max_cols=20
    )
    public_private_plot.save()  # TODO:
