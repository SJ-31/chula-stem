#!/usr/bin/env ipython

import io
import re
import subprocess
from collections.abc import Sequence
from functools import reduce
from pathlib import Path
from subprocess import run
from tempfile import NamedTemporaryFile
from typing import Literal

import anndata as ad
import awkward as ak
import matplotlib
import matplotlib.pyplot as plt
import mudata as md
import numpy as np
import pandas as pd
import plotnine as gg
import polars as pl
import polars.selectors as cs
import scirpy as ir
import seaborn as sns
import yaml
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from plotnine.ggplot import ggplot
from pyhere import here
from scipy import stats

with open(here("analyses", "pdac_tcr", "env.yaml"), "r") as f:
    config = yaml.safe_load(f)

source("/home/shannc/Bio_SDD/chula-stem/snakemake/tcr/prepare_data.py")
source("/home/shannc/Bio_SDD/chula-stem/snakemake/tcr/alignment.py")
source("/home/shannc/Bio_SDD/chula-stem/snakemake/tcr/plotting.py")
source("/home/shannc/Bio_SDD/chula-stem/snakemake/tcr/analyses.py")
# %%


run_name = "TCR2"

curdir = here("analyses", "pdac_tcr")

config["bd_results"][run_name] = here("output", "PDAC_TCR_2025-09-24/TCR2")

# %%


tag_mapping: dict = {f"SampleTag{i}_hs": v for i, v in config["sample_tags"].items()}

ct_col = "Cell_Type_Experimental"

# %%

mdata = md.read_h5mu(
    "/home/shannc/Bio_SDD/chula-stem/analyses/output/pdac_tcr/no_date/combined.h5mu"
)
mixcr_file = here("analyses", "output", "pdac_tcr", "testzone", "mixcr_PDAC82_airr.tsv")
mixcr2_file = here(
    "analyses", "output", "pdac_tcr", "testzone", "mixcr_PDAC82_tso_airr.tsv"
)


mixcr = ir.io.read_airr(mixcr_file)
ir.pp.index_chains(mixcr)
ir.tl.chain_qc(mixcr)
ir.tl.define_clonotypes(mixcr)

mixcr2 = ir.io.read_airr(mixcr2_file)
ir.pp.index_chains(mixcr2)
ir.tl.chain_qc(mixcr2)
ir.tl.define_clonotypes(mixcr2)

airr = mdata["airr"]
sample = "PDAC82"


ori = airr[airr.obs["Sample_Name"] == sample, :]

(~np.any(ir.get.airr(ori, "cdr3").iloc[:, :2].isna().values, axis=1)).sum()

(~np.any(ir.get.airr(ori, "cdr3").iloc[:, :2].isna().values, axis=1)).sum()


def compare_top_clones(x, y):
    compare_cols: list = ["v_call", "j_call", "c_call", "junction", "d_call"]
    tops: dict = {}
    for cur, adata in zip(["x", "y"], [x, y]):
        top = (
            adata.obs.loc[:, ["clone_id_size", "clone_id"]]
            .drop_duplicates()
            .sort_values("clone_id_size")
            .iloc[-1, :]
        )
        print(top)
        airr_data = ir.get.airr(adata, compare_cols)
        top_airr = airr_data.loc[top.name, :]
        tops[cur] = top_airr
    return tops


mixcr_comp = compare_top_clones(mixcr, ori)


airrs = {
    s: mdata[mdata.obs["Sample_Name"] == s, :]
    for s in mdata.obs["Sample_Name"].unique()
}

_, filtered = maybe_filter_by_rank(ori, config["mixtcrpred"])


cols = ["sequence", "cdr3", "cdr1", "cdr2", "fwr1", "fwr2", "fwr3", "fwr4"] + [
    s
    for seq in [
        [f"{k}_sequence_start", f"{k}_sequence_end", f"{k}_call"]
        for k in ("v", "d", "j")
    ]
    for s in seq
]

# Gonna need to align the cdrs and the framework regions

seqs = pl.from_pandas(ir.get.airr(filtered, ["sequence_id"] + cols + ["c_call"]))
vj_1 = seqs.select(cs.starts_with("VJ_1"))

align_test = Path(
    "/home/shannc/Bio_SDD/chula-stem/analyses/output/pdac_tcr/testzone/alignments"
)
inter_test = Path(
    "/home/shannc/Bio_SDD/chula-stem/analyses/output/pdac_tcr/testzone/alignments_inter"
)

align_vdj(vj_1, "VJ_1", align_test, id_col="VJ_1_sequence_id")

# seq_lengths =
#

vj_1["VJ_1_j_sequence_start"].mean()
vj_1["VJ_1_sequence"].str.len_chars().mean()

vj_1.with_columns(pl.col("VJ_1_sequence").str.slice(-pl.col("VJ_1_j_sequence_start")))

fig = plot_vdj(
    "10825037_TRA_1",
    "VJ_1",
    vj_1,
    align_test / "10825037_TRA_1.fasta",
    wrap_length=400,
    id_col="VJ_1_sequence_id",
    title_spec=[("", "VJ_1_sequence_id")],
    color_scheme="Clustal",
)
fig.set_figwidth(40)
fig.set_dpi(500)
fig.savefig(align_test.parent.joinpath("fig.png"))


# %%

wanted = ir.get.airr(mdata["airr"], "v_call")["VDJ_1_v_call"].unique()[:50]
leaders = get_leaders(wanted)

# %%
# * v gene usage

airr = mdata["airr"]
dfs = [
    ir.get.airr(airr[airr.obs["Sample_Name"] == s, :], "v_call")
    .loc[:, ["VJ_1_v_call", "VDJ_1_v_call"]]
    .melt()["value"]
    .value_counts()
    .rename(s)
    .reset_index()
    for s in airr.obs["Sample_Name"].unique()
]
# Want to determine if the v gene distribution is the same between samples
# So chi-squared, null hypothesis is that the distribution of each v gene is independent
# of sample
v_joined = (
    reduce(lambda x, y: x.merge(y, on="value", how="left"), dfs)
    .fillna(0)
    .set_index("value")
)
transposed = v_joined.transpose()
test = stats.chisquare(f_obs=transposed.iloc[:3, :])  # TODO: iloc isn't needed if you
# don't include the nonsense samples
v_joined["p_value"] = test.pvalue
# v_joined["p_adj"] =

mixtest = pd.read_csv(
    "/home/shannc/Bio_SDD/stem_synology/chula_mount/shannc/repos/MixTCRpred/test/test.csv"
)

val = format_tcrdb(
    "foo", Path("/home/shannc/Bio_SDD/chula-stem/analyses/data/tcrdb/immunoSEQ72")
)

# * Stitchr output


def prepare_for_thimble(
    data: md.MuData | ad.AnnData,
    receptor_class: Literal["AB", "GD"],
    multiple: bool = True,
    airr_mod: str = ["airr"],
) -> tuple[pl.DataFrame, pl.DataFrame]:
    data = data[airr_mod] if isinstance(data, md.MuData) else data
    target = "TRA+TRB" if receptor_class == "AB" else "TRG+TRD"
    obs: pl.DataFrame = pl.from_pandas(
        data.obs.reset_index().loc[:, ["chain_pairing", "receptor_subtype"]]
    )
    cdr3 = cdr
    airr_data = ir.get.airr(data, ["v_call", "j_call", "cdr3"])


# will

mdata["airr"].obsm["airr"].type.show()


# * Rarefaction curve


# %%

# * Clonotype network analysis
#
# ir.tl.clonotype_network(mdata, min_cells=10)
# ir.tl.repertoire_overlap(mdata, "Sample_Name")
# ir.tl.alpha_diversity(mdata, "Sample_Name", metric="pielou_e")
# ir.tl.alpha_diversity(mdata, "Sample_Name", metric="simpson")


# * Query reference dbs

vjdb = ir.datasets.vdjdb(cached=True, cache_path=here("analyses", "data", "vdjdb.h5ad"))
vjdb.obs_names = vjdb.obs_names + "_vjdb"
iedb = ir.datasets.iedb(cached=True, cache_path=here("analyses", "data", "iedb.h5ad"))
vjdb.obs = vjdb.obs.rename(columns=lambda x: f"vjdb_{x}")
mcpas_df = pl.read_csv(
    "/home/shannc/Bio_SDD/chula-stem/analyses/data/McPAS-TCR.csv",
    infer_schema_length=None,
)

cedar = pl.read_csv("/home/shannc/Downloads/tcr_full_v3.csv")

sources = (
    Path("/home/shannc/Downloads/tcr_full_v3.csv")
    .read_text()
    .splitlines()[0]
    .replace('"', "")
    .replace(" ", "_")
    .split(",")
)

cedar = cedar.rename()
name_mapping = {
    o: f"{s}_{n.replace(" ", "_")}"
    for o, n, s in zip(cedar.columns, cedar.row(0), sources)
}


# TODO: parse cedar database into airr format
def format_cedar(path) -> ad.AnnData:
    tmp = pl.read_csv(path, infer_schema_length=None)
    sources = (
        Path(path)
        .read_text()
        .splitlines()[0]
        .replace('"', "")
        .replace(" ", "_")
        .split(",")
    )
    name_mapping = {
        o: f"{s}_{n.replace(" ", "_")}"
        for o, n, s in zip(tmp.columns, tmp.row(0), sources)
    }
    tmp = (
        tmp.rename(name_mapping)
        .slice(1)
        .with_columns(
            pl.col("Receptor_Type").replace(
                {"alphabeta": "TRA+TRB", "gammadelta": "TRG+TRD"}
            )
        )
    )
    cedar_cells: list = []
    col_map: dict = {
        "v_call": "Curated_V_Gene",
        "j_call": "Curated_J_Gene",
        "d_call": "Curated_D_Gene",
        "sequence_aa": "Protein_Sequence",
        "sequence": "Nucleotide_Sequence",
        "junction_aa": "CDR3_Calculated",
    }
    for i, row in enumerate(tmp.iter_rows(named=True)):
        cell = ir.io.AirrCell(cell_id=i)
        chain_1 = ir.io.AirrCell.empty_chain_dict()
        chain_2 = ir.io.AirrCell.empty_chain_dict()

        chain_1_vals = {k: row[f"Chain_1_{v}"] for k, v in col_map.items()}
        chain_1.update({"productive": True, "locus": "TRA", "consensus_count": 0})
        chain_1.update(chain_1_vals)

        chain_2_vals = {k: row[f"Chain_2_{v}"] for k, v in col_map.items()}
        chain_2.update({"productive": True, "locus": "TRB", "consensus_count": 0})
        chain_2.update(chain_2_vals)

        cell.add_chain(chain_1)
        cell.add_chain(chain_2)
        cedar_cells.append(cell)
    cedar = ir.io.from_airr_cells(cedar_cells)
    # cedar.obs = tmp.drop(
    #     [
    #         "TRBV",
    #         "TRBJ",
    #         "TRBD",
    #         "TRAV",
    #         "TRAJ",
    #         "CDR3.alpha.nt",
    #         "CDR3.alpha.aa",
    #         "CDR3.beta.nt",
    #         "CDR3.beta.aa",
    #     ]
    # ).to_pandas()
    cedar.uns["DB"] = {"name": "McPAS", "date": None}
    ir.pp.index_chains(cedar)
    return cedar


# mcpas = format_mcpas("/home/shannc/Bio_SDD/chula-stem/analyses/data/McPAS-TCR.csv")

# pl.from_pandas(ir.get.airr(vjdb, "junction_aa")).with_columns(
#     pl.all().map_elements(lambda x: len(x) if x else 0, return_dtype=pl.Int64)
# ).max()

# pl.from_pandas(ir.get.airr(mcpas, "junction_aa")).with_columns(
#     pl.all().map_elements(lambda x: len(x) if x else 0, return_dtype=pl.Int64)
# ).max()
