#!/usr/bin/env ipython
from collections.abc import Sequence
from functools import reduce
from pathlib import Path
from typing import Literal

import anndata as ad
import awkward as ak
import matplotlib.pyplot as plt
import polars as pl
import scirpy as ir
import seaborn as sns
from matplotlib.figure import Figure
from snakemake.script import snakemake as smk

# * Helper functions


def format_mcpas(path) -> ad.AnnData:
    tmp = pl.read_csv(path, infer_schema_length=None)
    for i, row in enumerate(tmp.iter_rows(named=True)):
        cell = ir.io.AirrCell(cell_id=i)
        alpha_chain = ir.io.AirrCell.empty_chain_dict()
        beta_chain = ir.io.AirrCell.empty_chain_dict()
        alpha_chain.update(
            {
                "productive": True,
                "locus": "TRA",
                "consensus_count": 0,
                "junction": row["CDR3.alpha.nt"],
                "junction_aa": row["CDR3.alpha.aa"],
                "v_call": row["TRAV"],
                "j_call": row["TRAJ"],
                "d_call": None,
            }
        )
        beta_chain.update(
            {
                "productive": True,
                "locus": "TRB",
                "consensus_count": 0,
                "junction": row["CDR3.beta.nt"],
                "junction_aa": row["CDR3.beta.aa"],
                "v_call": row["TRBV"],
                "j_call": row["TRBJ"],
                "d_call": row["TRBD"],
            }
        )
        cell.add_chain(alpha_chain)
        cell.add_chain(beta_chain)
        mcpas_cells.append(cell)
    mcpas = ir.io.from_airr_cells(mcpas_cells)
    mcpas.obs = tmp.drop(
        [
            "TRBV",
            "TRBJ",
            "TRBD",
            "TRAV",
            "TRAJ",
            "CDR3.alpha.nt",
            "CDR3.alpha.aa",
            "CDR3.beta.nt",
            "CDR3.beta.aa",
        ]
    ).to_pandas()
    ir.pp.index_chains(mcpas)
    return mcpas


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


# * Load data from Rhapsody and save to a single h5ad object
if smk.rule == "load_runs":
    tag_mapping: dict = {
        f"SampleTag{i}_hs": v for i, v in smk.config["sample_tags"].items()
    }
    to_plot = smk.config["to_plot"]
    ct_col = "cell_type_experimental"

    adatas = []
    for run_name, dir in smk.config["bd_results"].items():
        run: Path = Path(dir)
        adata = ir.io.read_airr(
            run.joinpath(f"{run_name}_VDJ_Dominant_Contigs_AIRR.tsv")
        )
        dfs = []
        # TODO: this should be optional if samples were not multiplexed
        sample_df: dict = pl.read_csv(
            run.joinpath(f"{run_name}_Sample_Tag_Calls.csv"), comment_prefix="#"
        ).with_columns(pl.col("Sample_Name").replace(tag_mapping))
        dfs.append(sample_df)

        for dimr in ["tSNE", "UMAP"]:
            dimr_file = run.joinpath(f"{run_name}_{dimr}_coordinates.csv")
            if dimr_file.exists():
                dfs.append(pl.read_csv(dimr_file))

        if len(dfs) > 1:
            obs: pl.DataFrame = reduce(
                lambda x, y: x.join(y, on="Cell_Index", how="left"), dfs
            )
        else:
            obs = sample_df
        obs = (
            obs.with_columns(pl.col("Cell_Index").cast(pl.String))
            .filter(pl.col("Cell_Index").is_in(adata.obs.index))
            .with_columns(
                pl.lit(run_name).alias("run"),
                pl.Series(ak.to_list(adata.obsm["airr"][ct_col][:, 0])).alias(ct_col),
            )
        )
        adata.obs = adata.obs.merge(
            obs.to_pandas(),
            left_index=True,
            right_on="Cell_Index",
            how="left",
        ).set_index("Cell_Index", drop=False)

        ir.pp.index_chains(adata)
        ir.tl.chain_qc(adata)

        # Plot on per-run basis
        tsne_fig = plot_obs(adata, "tSNE", hues=to_plot)
        umap_fig = plot_obs(adata, "UMAP", hues=to_plot)

        adatas.append(adata)

    combined: ad.AnnData = ad.concat(adatas)

    # ** Call clonotypes
    # TODO: verify that you can do this with multiple samples
    calling_config = smk.config["clonotype_calling"]

    ir.pp.ir_dist(combined, **calling_config["ir_dist"])
    ir.tl.define_clonotypes(combined, **calling_config["define_clonotypes"])
    ir.tl.define_clonotype_clusters(
        combined, key_added="cc", **calling_config["define_clonotype_clusters"]
    )

    ir.tl.clonal_expansion(
        combined, expanded_in="Sample_Name", **smk.config["clonal_expansion"]
    )

    combined.obs_names_make_unique()
elif smk.rule == "prepare_references":
    mcpas_cells = []
    for k, v in smk.params["refs_to_get"].items():
        if k == "mcpas":
            adata = format_mcpas(path=v.parent.joinpath(f"{v.stem}.csv"))
            adata.write_h5ad(v)
        elif k == "vdjdb":
            _ = ir.datasets.vdjdb(cached=True, cache_path=v)
        elif k == "iedb":
            _ = ir.datasets.iedb(cached=True, cache_path=v)
