#!/usr/bin/env ipython

import anndata as ad
import pandas as pd
import plotnine as gg

try:
    from snakemake.script import snakemake as smk
except ImportError:
    smk = {}


def plot_group_abundance(
    adata: ad.AnnData,
    x: str,
    fill: str = "clonal_expansion",
    normalize: bool = False,
    max_cols: int = 50,
):
    "Custom plot because the scirpy version won't show legends for some reason"
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
        + gg.ylab("Number of Cells")
        + gg.xlab(x)
    )
    return plot


# * Generate plots
if smk.rule == "make_plots":
    adata = ad.read_h5ad(smk.input[0])
    ct_col = smk.config["ct_col"]
    airrs = {
        s: adata[adata.obs["Sample_Name"] == s, :]
        for s in adata.obs["Sample_Name"].unique()
        if s not in {"Undetermined", "Multiplet"}
    }
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
