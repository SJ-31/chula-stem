#!/usr/bin/env ipython
from collections.abc import Sequence
from pathlib import Path

import anndata as ad
import joblib
import mudata as md
import pandas as pd
import scirpy as ir

try:
    from snakemake.script import snakemake as smk
except ImportError:
    smk = type("snakemake", (), {"rule": None})


def query_routine(
    data: md.MuData | ad.AnnData,
    reference: md.MuData | ad.AnnData,
    include_query_cols: Sequence | None = None,
    include_ref_cols: Sequence | None = None,
    **kws,
) -> pd.DataFrame:
    """Wrapper for scirpy's ir query pipeline
    Results are a dataframe from ir_query_annotate_df

    Note: This modifies `data` inplace

    Parameters
    ----------
    kws : Keyword arguments passed to pp.ir_dist and tl.ir_query
    """
    dist_kws: dict = {k: v for k, v in kws.items() if k in {"metric", "sequence"}}
    cutoff = kws.pop("cutoff", None)
    match_columns = kws.pop("match_columns", None)
    ir.pp.ir_dist(data, reference, cutoff=cutoff, **dist_kws)
    ir.tl.ir_query(
        data,
        reference,
        match_columns=match_columns,
        **kws,
    )
    match_df = ir.tl.ir_query_annotate_df(
        data,
        reference,
        include_query_cols=include_query_cols,
        include_ref_cols=include_ref_cols,
        **dist_kws,
    )
    return match_df


# * Entry point
if smk.rule == "scirpy_query":
    mdata = md.read_h5mu(smk.input["mdata"])
    for db_name, path in smk.params["references"].items():
        reference = md.read_h5mu(path)
        outfile = Path(smk.params["outdir"]) / f"{db_name}_query"
        result: pd.DataFrame = query_routine(
            mdata,
            reference=reference,
            include_query_cols=["Sample_Name", "airr:clone_id", "airr:clone_id_size"],
            **smk.config["query_reference"]["scirpy"],
        )
        result.to_csv(outfile.with_suffix(".csv"), index=False)
        joblib.dump(mdata.uns, outfile.with_suffix(".pkl"))
