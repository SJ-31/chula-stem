#!/usr/bin/env ipython

from pathlib import Path
from subprocess import run

import numpy as np
import pandas as pd
import scirpy as ir

try:
    from snakemake.script import snakemake as smk
except ImportError:
    smk = type("snakemake", (), {"rule": None, "config": {}})

from analyses import SCOL, get_airr, maybe_filter_by_rank

RENAME: dict = {
    "VDJ_1_junction_aa": "cdr3_TRB",
    "VJ_1_junction_aa": "cdr3_TRA",
    "VDJ_1_v_call": "TRBV",
    "VDJ_1_j_call": "TRBJ",
    "VJ_1_v_call": "TRAV",
    "VJ_1_j_call": "TRAJ",
}

# * Prepare input data

airr = get_airr(filter_samples=True)
_, airr = maybe_filter_by_rank(airr, smk.config["mixtcrpred"])
mixtcr_config: dict = smk.config["mixtcrpred"]
with ir.get.airr_context(airr, ["junction_aa", "v_call", "j_call"]):
    input_obs: pd.DataFrame = airr.obs.rename(RENAME, axis=1)
length_mask = (input_obs["cdr3_TRA"].str.len() > 20) | (
    input_obs["cdr3_TRB"].str.len() > 20
)
input_obs = input_obs.loc[~length_mask, :]
if mixtcr_config.get("ignore_incomplete", False):
    input_obs = input_obs.loc[np.any(input_obs.loc[:, RENAME.values()].isna()), :]
input_obs = input_obs.drop_duplicates([SCOL] + list(RENAME.values()))

# * Run MixTCR


if (env := mixtcr_config.get("env")) and (conda_dir := smk.config.get("conda")):
    conda_flags = f"--conda_env {env} --conda_dir {conda_dir}"
else:
    conda_flags = ""
rundir = mixtcr_config["dir"]
if not Path(rundir).exists():
    raise ValueError("Path to MixTCRpred repo doesn't exist")
if not (meta_path := Path(rundir) / "pretrained_models" / "info_models.csv").exists():
    raise ValueError("MixTCRpred pretrained_models/info_models.csv file missing")
metadata: pd.DataFrame = pd.read_csv(meta_path)
input_file = Path(f"{smk.params["outdir"]}/input.csv")

dfs = []
input_obs.loc[:, RENAME.values()].reset_index().to_csv(input_file, index=False)
for model in mixtcr_config["models"]:
    if len(model.split("_")) != 2:
        raise ValueError("Invalid model specification for MixTCRpred!")
    outfile = f"{smk.params["outdir"]}/{model}.csv"
    command = (
        f"./mixtcr_wrapper.sh {rundir} {model} {input_file} {outfile} {conda_flags}"
    )
    run(command, shell=True)
    df = pd.read_csv(outfile, comment="#")
    df = (
        df.merge(input_obs.reset_index(), on=input_obs.index.name)
        .loc[
            :,
            [input_obs.index.name, SCOL, "score", "perc_rank"]
            + smk.config["obs_cols_include"],
        ]
        .assign(MixTCRpred_model_name=model)
        .merge(metadata, on="MixTCRpred_model_name")
    )
    dfs.append(df)
pd.concat(dfs).to_csv(smk.output[0], index=False)
