#!/usr/bin/env ipython
import io
import re
from collections.abc import Sequence
from pathlib import Path
from subprocess import run
from tempfile import NamedTemporaryFile

import anndata as ad
import joblib
import mudata as md
import numpy as np
import pandas as pd
import scirpy as ir

md.set_options(pull_on_update=False)

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


def get_airr(input_key: str | int = 0, filter_samples: bool = False) -> ad.AnnData:
    airr: ad.AnnData = md.read_h5mu(smk.input[0])["airr"]
    if filter_samples:
        airr = airr[~airr.obs[SCOL].isin(smk.config["ignore_samples"]), :]
    return airr


def maybe_filter_by_rank(adata: ad.AnnData, config: dict) -> tuple[bool, ad.AnnData]:
    if (rank_key := config.get("rank_key")) and (top := config.get("top")):
        return True, adata[adata.obs[rank_key] <= top, :]
    return False, adata


def fasta_from_airr(
    adata: ad.AnnData,
    min_length: int = 10,
    chains: tuple = ("VJ_1", "VDJ_1"),
    cols_add: tuple = ("rank", "sample"),
    return_aa: bool = False,
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
    return_aa : bool
        Extracted sequences are amino acids

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
    ignored.update({c: [] for c in cols_add})
    seq_str = "sequence" if not return_aa else "sequence_aa"
    cols = [
        "sequence_id",
        seq_str,
        "complete_vdj",
        "v_call",
        "c_call",
        "d_call",
        "j_call",
    ]
    for (id, obs_row), (_, airr_row) in zip(
        adata.obs.iterrows(), ir.get.airr(adata, cols).iterrows()
    ):
        subtype = obs_row["receptor_subtype"]
        clone_id = obs_row["clone_id"]
        is_ambiguous: bool = subtype == "ambiguous"
        no_ir: bool = subtype == "no IR"
        for chain in chains:
            sequence = airr_row[f"{chain}_{seq_str}"]
            seqid = airr_row[f"{chain}_sequence_id"]
            call_keys = ["v" "d", "j", "c"] if "D" in chains else ["v", "j", "c"]
            is_complete: bool = airr_row[f"{chain}_complete_vdj"]
            has_sequence: bool = sequence and len(sequence) > min_length
            ignore_if = [not is_complete, is_ambiguous, not has_sequence, no_ir]
            if any(ignore_if):
                ignored["clone_id"].append(clone_id)
                ignored["chain"].append(chain)
                ignored["receptor_subtype"].append(subtype)
                ignored["cell_id"].append(id)
                for c in cols_add:
                    ignored[c].append(obs_row[c])
                for reason, val in zip(
                    ["incomplete", "is_ambiguous", "no_sequence", "no_ir"], ignore_if
                ):
                    if val:
                        ignored["reason"].append(reason)
                        break
                continue
            calls = " ".join(
                [f"{k}_call={airr_row[f'{chain}_{k}_call']}" for k in call_keys]
            )
            chain_name = CHAIN_NAMING[re.sub("_[0-9]+.*", "", chain)].get(subtype, "na")
            header = f"{seqid} {calls} clone_id={clone_id} chain={chain_name} subtype={subtype}"
            for col in cols_add:
                header = f"{header} {col}={obs_row[col]}"
            sequences.append(f">{header}\n{sequence}")
    return "\n".join(sequences), pd.DataFrame(ignored)


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
    cols = ["v_call", "d_call", "j_call", "junction_aa"]
    query_seqs = ir.get.airr(data, cols).rename(lambda x: f"query_{x}", axis=1)
    ref_seqs = ir.get.airr(reference, cols).rename(lambda x: f"ref_{x}", axis=1)
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
    ).drop_duplicates()
    match_df = match_df.merge(query_seqs, left_index=True, right_index=True, how="left")
    match_df = match_df.merge(ref_seqs, left_index=True, right_index=True, how="left")
    return match_df


# * Wrapper functions


def scirpy_query(airr: ad.AnnData, query_cols: list):
    for db_name, path in smk.params["references"].items():
        reference = ad.read_h5ad(path)
        outfile = Path(smk.params["outdir"]) / f"{db_name}_query"
        try:
            result: pd.DataFrame = query_routine(
                airr,
                reference=reference,
                include_query_cols=query_cols,
                **config["scirpy"],
            )
        except KeyError:  # BUG: something is wrong with the McPAS database
            result = pd.DataFrame()
        result.to_csv(outfile.with_suffix(".csv"), index=False)
        joblib.dump(airr.uns, outfile.with_suffix(".pkl"))


def tcrmatch_wrapper(
    database: str,
    adata: ad.AnnData | None = None,
    sequences: Sequence | Path | str | None = None,
    query_cols: tuple = ("clone_id", "clone_id_size"),
    memory: int = 4,
    threads: int = 1,
    threshold: float = 0.97,
    cdr3_col: str = "VDJ_1_cdr3_aa",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Wrapper function for tcrmatch

    Parameters
    ----------
    adata : ad.AnnData | None
        adata object containing air data
    sequences : Sequence | Path | str | None
        path to cdr3 beta sequences to analyze
    """
    seqs: Sequence
    invalid: pd.DataFrame = pd.DataFrame()
    tmp: pd.DataFrame = pd.DataFrame()
    query_cols = (
        ["clone_id"] + list(query_cols) if "clone_id" not in query_cols else query_cols
    )
    if adata is not None:
        wrong_subtype = adata[adata.obs["receptor_subtype"] != "TRA+TRB", :]
        valid = adata[adata.obs["receptor_subtype"] == "TRA+TRB", :]
        cdr3s = ir.get.airr(valid, "cdr3_aa")
        no_cdr3beta = valid[~cdr3s[cdr3_col].notnull(), :]
        invalid = pd.concat(
            [
                df.obs.loc[:, query_cols].reset_index().assign(reason=r)
                for df, r in zip(
                    [wrong_subtype, no_cdr3beta], ["wrong_subtype", "no_cdr3_beta"]
                )
            ]
        )
        mask = ~cdr3s.index.isin(no_cdr3beta.obs.index)
        tmp = (
            valid.obs.loc[mask, query_cols]
            .reset_index()
            .merge(cdr3s.loc[mask, [cdr3_col]].reset_index(), on="index")
        )
        tmp = tmp.loc[~tmp.duplicated([cdr3_col, "clone_id"]), :]
        seqs = tmp[cdr3_col]
    elif isinstance(sequences, Path):
        seqs = sequences.read_text().splitlines()
    elif isinstance(sequences, str):
        with open(sequences, "r") as f:
            seqs = f.read().splitlines()
    else:
        raise ValueError("Either `adata` or `sequences` must be provided!")
    database = str(database) if isinstance(database, Path) else database
    args = ["-d", database, "-t", str(threads), "-m", str(memory), "-s", str(threshold)]
    with NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write("\n".join(seqs))
    command = f"tcrmatch -i {f.name} {' '.join(args)}"
    stdout = run(command, shell=True, capture_output=True).stdout
    Path(f.name).unlink()
    if stdout is None:
        return pd.DataFrame(), invalid
    match_result: pd.DataFrame = pd.read_csv(io.StringIO(stdout.decode()), sep="\t")
    if match_result.empty:
        return pd.DataFrame(), invalid
    if adata is not None:
        match_result = (
            tmp.merge(
                match_result,
                left_on=cdr3_col,
                right_on="trimmed_input_sequence",
                how="right",
            )
            .rename({"index": "cell_id"}, axis=1)
            .drop(cdr3_col, axis=1)
        )
    return match_result, invalid


# * Entry point
# * Querying
if smk.rule in {"scirpy_query", "tcrmatch"}:
    airr = get_airr(0, filter_samples=True)
    config: dict = smk.config["query_reference"]
    query_cols = [SCOL, "clone_id", "clone_id_size"]
    filtered, airr = maybe_filter_by_rank(airr, config)
    if filtered:
        query_cols.append(config["rank_key"])
    # ** Scirpy
    if smk.rule == "scirpy_query":
        scirpy_query(airr, query_cols)
    # ** tcrmatch
    elif smk.rule == "tcrmatch":
        tcrmatch_config: dict = config["tcrmatch"]
        for i, (name, db) in enumerate(tcrmatch_config["databases"].items()):
            result, invalid = tcrmatch_wrapper(
                adata=airr,
                database=db,
                query_cols=query_cols,
                memory=tcrmatch_config.get("m", 4),
                threads=tcrmatch_config.get("threads", 1),
                threshold=tcrmatch_config.get("threshold", 0.97),
            )
            result.to_csv(f"{smk.params['outdir']}/tcrmatch_{name}.csv", index=False)
            if i == 0:
                invalid.to_csv(
                    f"{smk.params['outdir']}/tcrmatch_invalid.csv", index=False
                )
elif smk.rule == "extract_sequences":
    airr: ad.AnnData = md.read_h5mu(smk.input[0])["airr"]
    samples = set(airr.obs[SCOL].unique()) - set(smk.config.get("ignore_samples", []))
    config = smk.config.get("get_sequences", {})
    top = config.get("top", 3)
    rank_col = config.get("rank_key", "clone_id_Sample_Name_rank")
    airr.obs = airr.obs.rename(
        {rank_col: "rank", SCOL: "sample", "clone_id_size": "clone_size"}, axis=1
    )
    all_ignored = []
    for return_aa, suffix in zip([False, True], ["", "-aa"]):
        for i, sample in enumerate(samples):
            mask = (airr.obs["sample"] == sample) & (airr.obs["rank"] <= top)
            cur = airr[mask, :]
            fasta, ignored = fasta_from_airr(
                cur,
                min_length=config.get("min_length", 2),
                chains=("VJ_1", "VDJ_1"),
                cols_add=("rank", "sample", "clone_size"),
                return_aa=return_aa,
            )
            outfile = Path(smk.params["outdir"]) / f"{sample}{suffix}.fasta"
            outfile.write_text(fasta)
            if i == 0:
                all_ignored.append(ignored)
    pd.concat(all_ignored).to_csv(smk.output["ignored"], index=False)
