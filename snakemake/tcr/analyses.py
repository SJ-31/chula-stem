#!/usr/bin/env ipython
import re
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
    smk = type("snakemake", (), {"rule": None, "config": {}})

CHAIN_NAMING = {
    "VJ": {"TRG+TRD": "gamma", "TRA+TRB": "alpha"},
    "VDJ": {"TRG+TRD": "delta", "TRA+TRB": "beta"},
}
# TODO: you can include the other receptor subtypes  (BCRs mainly) if needed

SCOL: str = smk.config.get("sample_col", "Sample_Name")


def fasta_from_airr(
    adata: ad.AnnData,
    min_length: int = 10,
    chains: tuple = ("VJ_1", "VDJ_1"),
    cols_add: tuple = ("rank", "sample"),
) -> tuple[str, pd.DataFrame]:
    """Extract full-length air sequences from IR data in `adata`

    Parameters
    ----------
    min_length : int
        Minimum length of air sequence to allow
    chains : tuple
        Tuple listing air chain types to get sequence data from.
        Is the first chain pair by default
    cols_add : tuple
        Columns in adata.obs to annotate the sequence fasta header with, in the format of
        ... <column_name>=<value>


    Returns
    -------
    A tuple with elements
    1. valid airr sequences in `adata`, in a FASTA-formatted string
    2. dataframe with metadata about sequences that could not be obtained
    """
    sequences = []
    ignored = {
        "clone_id": [],
        "cell_id": [],
        "reason": [],
        "chain": [],
        "receptor_subtype": [],
    }
    cols = ["sequence_id", "sequence", "complete_vdj"]
    for (id, obs_row), (_, airr_row) in zip(
        adata.obs.iterrows(), ir.get.airr(adata, cols).iterrows()
    ):
        subtype = obs_row["receptor_subtype"]
        clone_id = obs_row["clone_id"]
        is_ambiguous: bool = subtype == "ambiguous"
        no_ir: bool = subtype == "no IR"
        for chain in chains:
            sequence = airr_row[f"{chain}_sequence"]
            seqid = airr_row[f"{chain}_sequence_id"]
            is_complete: bool = airr_row[f"{chain}_complete_vdj"]
            has_sequence: bool = sequence and len(sequence) > min_length
            ignore_if = [not is_complete, is_ambiguous, not has_sequence, no_ir]
            if any(ignore_if):
                ignored["clone_id"].append(clone_id)
                ignored["chain"].append(chain)
                ignored["receptor_subtype"].append(subtype)
                ignored["cell_id"].append(id)
                for reason, val in zip(
                    ["incomplete", "is_ambiguous", "no_sequence", "no_ir"], ignore_if
                ):
                    if val:
                        ignored["reason"].append(reason)
                        break
                continue
            chain_name = CHAIN_NAMING[re.sub("_[0-9]+.*", "", chain)].get(subtype, "na")
            header = f"{seqid} clone_id={clone_id} chain={chain_name} subtype={subtype}"
            for col in cols_add:
                header = f"{header} {col}={obs_row[col]}"
            sequences.append(f">{header}\n{sequence}")
    return "\n".join(sequences), pd.DataFrame(ignored, strict=False)


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
            include_query_cols=[SCOL, "airr:clone_id", "airr:clone_id_size"],
            **smk.config["query_reference"]["scirpy"],
        )
        result.to_csv(outfile.with_suffix(".csv"), index=False)
        joblib.dump(mdata.uns, outfile.with_suffix(".pkl"))
elif smk.rule == "extract_sequences":
    airr: ad.AnnData = md.read_h5mu(smk.input[0])["airr"]
    samples = set(airr.obs[SCOL].unique()) - set(smk.config.get("ignore_samples", []))
    config = smk.config("get_sequences", {})
    top = config.get("top", 3)
    rank_col = config.get("rank_key", "clone_id_Sample_Name_rank")
    airr.obs = airr.obs.rename({rank_col: "rank", SCOL: "sample"}, axis=1)
    all_ignored = []
    for sample in samples:
        mask = (airr.obs["sample"] == sample) & (airr.obs["rank"] <= top)
        cur = airr[mask, :]
        fasta, ignored = fasta_from_airr(
            cur,
            min_length=config.get("min_length", 2),
            chains=("VJ_1", "VDJ_1"),
            cols_add=("rank", "sample"),
        )
        all_ignored.append(ignored.assign(**{SCOL: sample}))
        outfile = Path(smk.params["outdir"]) / f"{sample}.fasta"
        outfile.write_text(fasta)
    pd.concat(all_ignored).to_csv(smk.output["ignored"], index=False)
