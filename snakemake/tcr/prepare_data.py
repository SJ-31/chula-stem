#!/usr/bin/env ipython
from pathlib import Path

import anndata as ad
import mudata as md
import pandas as pd
import polars as pl
import scirpy as ir

try:
    from snakemake.script import snakemake as smk
except ImportError:
    smk = type("snakemake", (), {"rule": None})
md.set_options(pull_on_update=False)


# * Helper functions
def alpha_richness(
    data: ad.AnnData | md.MuData,
    groupby: str,
    airr_mod: str = "airr",
    target_col: str = "clone_id",
    normalize: bool = True,
) -> None:
    if isinstance(data, md.MuData):
        data = data[airr_mod]
    richness_map = data.obs.groupby(groupby, observed=True).agg(
        {target_col: pd.Series.nunique}
    )[target_col]
    if normalize:
        richness_map = richness_map / richness_map.sum()
    data.obs = data.obs.assign(
        richness_clone_id=data.obs[groupby].cat.rename_categories(
            richness_map.to_dict()
        )
    )


def mark_public_clones(
    adata: ad.AnnData, clone_col: str = "clone_id", sample_col: str = "Sample_Name"
) -> None:
    public_clones = (
        pl.from_pandas(adata.obs)
        .filter(
            ~pl.col(sample_col).is_in(["Multiplet", "Undetermined"])
            & pl.col(clone_col).is_not_null()
        )
        .group_by(clone_col)
        .agg(pl.col(sample_col).unique().len().alias("count"))
        .filter(pl.col("count") > 1)["clone_id"]
        .to_list()
    )
    adata.obs.loc[:, "public_clonotype"] = adata.obs[clone_col].isin(public_clones)


def format_mcpas(path) -> ad.AnnData:
    tmp = pl.read_csv(path, infer_schema_length=None)
    mcpas_cells: list = []
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


def load_rhapsody_run(path: Path, run_name: str, tag_mapping: dict) -> md.MuData:
    airr: ad.AnnData = ir.io.read_airr(
        path.joinpath(f"{run_name}_VDJ_Dominant_Contigs_AIRR.tsv")
    )
    rna: ad.AnnData = md.read_h5mu(path.joinpath(f"{run_name}.cellismo")).mod["rna"]
    rna.obs = rna.obs.assign(
        Sample_Name=rna.obs["Sample_Name"].cat.rename_categories(tag_mapping)
    )
    mdata: md.MuData = md.MuData({"rna": rna, "airr": airr})
    mdata.pull_obs(mods="rna", prefix_unique=False, drop=True)
    mdata.push_obs("Sample_Name", mods="airr")
    return mdata


# * Load data from Rhapsody and save to a single mudata object
def load_runs():
    tag_mapping: dict = {
        f"SampleTag{i}_hs": v for i, v in smk.config["sample_tags"].items()
    }
    scol = "Sample_Name"

    mdatas = []
    for run_name, dir in smk.params["runs"].items():
        run: Path = Path(dir)
        mdata = load_rhapsody_run(run, run_name, tag_mapping=tag_mapping)
        mdata.obs_names = [f"{x}_{run_name}" for x in mdata.obs_names]
        ir.pp.index_chains(mdata)
        ir.tl.chain_qc(mdata)
        mdatas.append(mdata)

    if len(mdatas) > 1:
        combined: md.MuData = md.concat([mdatas])
    else:
        combined = mdatas[0]

    # ** Call clonotypes
    # TODO: verify that you can do this with multiple samples
    calling_config = smk.config["clonotype_calling"]

    ir.tl.define_clonotypes(combined)
    shared_kws = {
        k: v
        for k, v in calling_config["define_clonotype_clusters"].items()
        if k in {"metric", "sequence"}
    }
    calling_config["ir_dist"].update(shared_kws)
    ir.pp.ir_dist(combined, **calling_config["ir_dist"])
    ir.tl.define_clonotype_clusters(
        combined, **calling_config["define_clonotype_clusters"]
    )
    ir.tl.clonal_expansion(combined, expanded_in=scol, **smk.config["clonal_expansion"])
    mark_public_clones(combined["airr"])

    # ** Diversity
    for alpha in smk.config["alpha_metrics"]:
        if alpha == "richness":
            alpha_richness(combined["airr"], groupby=scol, target_col="clone_id")
        else:
            ir.tl.alpha_diversity(combined, groupby=scol, metric=alpha)

    combined.write_h5mu(smk.output[0])


# * Prepare reference files
def prepare_references():
    for k, v in smk.params["refs_to_get"].items():
        path = Path(v)
        as_adata = path.with_suffix(".h5ad")
        if k == "mcpas":
            if as_adata.exists():
                adata = ad.read_h5ad(as_adata)
            else:
                adata = format_mcpas(path=path.with_suffix(".csv"))
                adata.write_h5ad(v)
        elif k == "vdjdb":
            adata = ir.datasets.vdjdb(cached=True, cache_path=as_adata)
        elif k == "iedb":
            adata = ir.datasets.iedb(cached=True, cache_path=as_adata)
        elif as_adata.exists():
            continue
        else:  # Read in an aribtary airr-formatted tsv file
            raise NotImplementedError("")


# * Custom annotation with CellAssign
def annotate_cells():
    import numpy as np
    import pandas as pd
    from scvi.external import CellAssign

    config = smk.config[smk.rule]["cellassign"]
    mdata: md.MuData = md.read_h5mu(smk.input[0])

    markers: pd.DataFrame = pd.read_csv(config["marker_file"]).set_index(
        config["symbol_col"]
    )

    lib_size = mdata["rna"].X.sum(1)
    mdata["rna"].obs["size_factor"] = lib_size / np.mean(lib_size)
    filtered = mdata["rna"][:, mdata["rna"].var.index.isin(markers.index)].copy()
    CellAssign.setup_anndata(filtered, size_factor_key="size_factor")

    model = CellAssign(filtered, markers)
    model.train(**config.get("train_kws", {}))

    predictions = model.predict()
    elbo_plot = model.history["elbo_validation"].plot()

    predictions.index = mdata.obs.index
    predictions.reset_index().to_csv(smk.output["predictions_raw"], index=False)
    elbo_plot.figure.savefig(smk.output["elbo_plot"])  # TODO:
    result = pl.DataFrame(
        {"index": filtered.obs.index, "prediction": predictions.idxmax(axis=1).values}
    )
    result.write_csv(smk.output["predictions"])


# * Entry point

if smk.rule == "load_runs":
    load_runs()
elif smk.rule == "annotate_cells":
    annotate_cells()
elif smk.rule == "prepare_references":
    prepare_references()
