#!/usr/bin/env ipython
from pathlib import Path

import anndata as ad
import mudata as md
import polars as pl
import scirpy as ir

try:
    from snakemake.script import snakemake as smk
except ImportError:
    smk = type("snakemake", (), {"rule": None})
md.set_options(pull_on_update=False)

# * Helper functions


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
    return mdata


# * Load data from Rhapsody and save to a single h5ad object
# TODO: need to include the expression data from cellismo
if smk.rule == "load_runs":
    tag_mapping: dict = {
        f"SampleTag{i}_hs": v for i, v in smk.config["sample_tags"].items()
    }
    to_plot = smk.config["to_plot"]
    ct_col = "cell_type_experimental"

    mdatas = []
    for run_name, dir in smk.config["bd_results"].items():
        run: Path = Path(dir)
        mdata = load_rhapsody_run(run, run_name, tag_mapping=tag_mapping)
        mdata.obs_names = [f"{x}_{run_name}" for x in mdata.obs_names]
        ir.pp.index_chains(mdata)
        ir.tl.chain_qc(mdata)
        mdatas.append(mdata)

    combined: md.MuData = md.concat([mdatas], index_unique=True)

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
    mark_public_clones(combined)
    # combined.write_h5ad(smk.output)


elif smk.rule == "prepare_references":
    mcpas_cells = []
    for k, v in smk.params["refs_to_get"].items():
        if k == "mcpas":
            mdata = format_mcpas(path=v.parent.joinpath(f"{v.stem}.csv"))
            mdata.write_h5ad(v)
        elif k == "vdjdb":
            _ = ir.datasets.vdjdb(cached=True, cache_path=v)
        elif k == "iedb":
            _ = ir.datasets.iedb(cached=True, cache_path=v)
