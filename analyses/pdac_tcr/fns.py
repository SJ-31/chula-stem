#!/usr/bin/env ipython

import anndata as ad
import mudata as md
import polars as pl
import scirpy as ir
import yaml
from pyhere import here

with open(here("analyses", "pdac_tcr", "env.yaml"), "r") as f:
    CONFIG = yaml.safe_load(f)


def get_airrs():
    mdata = md.read_h5mu(
        "/home/shannc/Bio_SDD/chula-stem/analyses/output/pdac_tcr/no_date/combined.h5mu"
    )
    mixcr_file = here(
        "analyses", "output", "pdac_tcr", "testzone", "mixcr_PDAC82_airr.tsv"
    )
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
    return {"ori": ori, "mixcr": mixcr}


def compare_null_counts(adatas: ad.AnnData, to_get, names):
    dfs = [pl.from_pandas(ir.get.airr(a, to_get)) for a, n in zip(adatas, names)]
    result = {"column": []}
    result.update(
        {f"{n}_{val}": [] for n in names for val in ["null_percent", "total"]}
    )
    for col in dfs[0].columns:
        result["column"].append(col)
        for d, n in zip(dfs, names):
            null_percent = d[col].is_null().sum() / d.shape[0] * 100
            result[f"{n}_total"].append(d.shape[0])
            result[f"{n}_null_percent"].append(null_percent)
    return pl.DataFrame(result)
