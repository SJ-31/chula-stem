#!/usr/bin/env ipython
import io
import re
from collections.abc import Sequence
from pathlib import Path
from subprocess import Popen, run
from tempfile import NamedTemporaryFile
from typing import Literal

import anndata as ad
import joblib
import mudata as md
import pandas as pd
import polars as pl
import polars.selectors as cs
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
    airr: ad.AnnData = md.read_h5mu(smk.input[input_key])["airr"]
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
    # BUG: [2025-10-10 Fri] this won't work for tcrdb data because it needs both VJ and VDJ chains...
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
        if "tcrdb" in db_name:
            continue  # query_routine doesn't work with orphan VDJ
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


def format_for_compairr(
    adata: ad.AnnData,
    repertoire: str,
    nucleotides: bool = False,
    cdr3: bool = False,
    repertoire_from_obs: bool = False,
    new_seqids: bool = False,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    cols = ["sequence_id", "v_call", "j_call", "junction", "junction_aa"]
    chains = ("VJ_1", "VDJ_1")
    iname: str = adata.obs.index.name
    df = pl.from_pandas(ir.get.airr(adata, cols), include_index=True)
    if cdr3:
        cols.extend(["cdr3", "cdr3_aa"])
    if new_seqids:
        exprs = [pl.col(iname).alias(f"{c}_sequence_id") for c in chains]
        df = df.with_columns(*exprs)
    suffix = "_aa" if not nucleotides else ""
    to_group_by = ["v_call", "j_call", f"junction{suffix}"]
    if cdr3:
        to_group_by.append(f"cdr3{suffix}")
    keep_first = (set(cols) - set(to_group_by)) | {"repertoire_id"}
    if not repertoire_from_obs:
        df = df.with_columns(pl.lit(repertoire).alias("repertoire_id"))
    else:
        df = df.join(
            pl.from_pandas(adata.obs.loc[:, [repertoire]], include_index=True), on=iname
        ).with_columns(pl.col(repertoire).alias("repertoire_id"))
    df = (
        pl.concat(
            [
                df.select(cs.starts_with(chain) | pl.col("repertoire_id"))
                .rename(lambda x: x.replace(f"{chain}_", ""))
                .with_columns(pl.col("sequence_id") + pl.lit(f"-{chain}"))
                for chain in chains
            ]
        )
        .sort("sequence_id")
        .group_by(to_group_by)
        .agg(pl.col(keep_first).first(), pl.len().alias("duplicate_count"))
        .with_columns(
            pl.all_horizontal(pl.col(to_group_by).is_not_null()).alias("all"),
        )
    )
    invalid = df.filter(pl.col("all").not_()).drop("all")
    return df.filter(pl.col("all")).drop("all"), invalid


def compairr_wrapper(
    query: ad.AnnData | md.MuData,
    reference: ad.AnnData | md.MuData,
    repertoire_col: str,
    ref_name: str = "",
    repertoire_from_obs: bool = True,
    nucleotides: bool = False,
    differences: int = 0,
    cdr3: bool = False,
    indels: bool = False,
    threads: int = 1,
    score: Literal["product", "ratio", "min", "max", "mean"] = "product",
    switches: Sequence = ("ignore-unknown",),
    include_ref_obs: bool = True,
    new_ref_seqids: bool = False,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Wrapper function for calling compairr in query mode (--existence switch)

    Parameters
    ----------
    param : argument
    switches : Sequence
        Additional switches to pass to compairr, WITHOUT the prefix dashes
        e.g. "ignore-counts", "distance"
    """
    valid, invalid = format_for_compairr(
        query,
        repertoire_col,
        nucleotides=nucleotides,
        cdr3=cdr3,
        repertoire_from_obs=repertoire_from_obs,
    )
    valid_r, _ = format_for_compairr(
        reference,
        ref_name,
        nucleotides=nucleotides,
        cdr3=cdr3,
        new_seqids=new_ref_seqids,
    )
    args = [
        "--existence",
        "--differences",
        differences,
        "--threads",
        threads,
        "--score",
        score,
    ]
    valid_switches = {"ignore-unknown", "ignore-genes", "ignore-counts", "indels"}
    for switch in switches:
        if switch in valid_switches:
            args.append(f"--{switch}")
        else:
            raise ValueError(f"Switch {switch} not supported by compairr")
    if cdr3:
        args.append("--cdr3")
    if nucleotides:
        args.append("--nucleotides")
    if indels:
        args.append("--indels")
    with NamedTemporaryFile("+w", suffix="tsv") as q:
        with NamedTemporaryFile("+w", suffix="tsv") as r:
            with NamedTemporaryFile("+w", suffix="tsv") as outfile:
                valid.write_csv(q.name, separator="\t")
                valid_r.write_csv(r.name, separator="\t")
                args.extend(["--outfile", outfile.name])
                args = [str(a) for a in args]
                with Popen(["compairr", q.name, r.name] + args) as proc:
                    _ = proc.communicate()
                    pl.read_csv(outfile, separator="\t")
    chain_exprs = [
        [
            pl.col(f"sequence_id_{c}").str.extract(r"-(VD?J_[12])").alias(f"chain_{c}"),
            pl.col(f"sequence_id_{c}").str.replace_many(["-VDJ_1", "-VJ_1"], ["", ""]),
        ]
        for c in ["query", "ref"]
    ]
    result = result.rename(
        lambda x: x.replace("_1", "_query").replace("_2", "_ref")
    ).with_columns(*chain_exprs[0], *chain_exprs[1])
    if include_ref_obs:
        id_col = "sequence_id_ref"
        ref_obs: pl.DataFrame = pl.from_pandas(reference.obs, include_index=True)
        result = result.join(ref_obs, left_on=id_col, right_on=reference.obs.index.name)
    return result, invalid


# * Querying
if smk.rule in {"scirpy_query", "tcrmatch"}:
    airr = get_airr(0, filter_samples=True)
    config: dict = smk.config["query_reference"]
    query_cols = [SCOL, "clone_id", "clone_id_size"]
    did_filter, airr = maybe_filter_by_rank(airr, config)
    if did_filter:
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
    # ** compairr
    elif smk.rule == "compairr":
        compairr_config: dict = config["compairr"]
        databases = smk.params["references"]
        for i, (name, db_path) in enumerate(databases.items()):
            result, invalid = compairr_wrapper(
                airr,
                ad.read_h5ad(db_path),
                repertoire_col=SCOL,
                ref_name=name,
                new_ref_seqids=True,
                **compairr_config,
            )
            result.write_csv(f"{smk.params['outdir']}/{name}.csv")
            if i == 0:
                invalid.write_csv(f"{smk.params['outdir']}/invalid.csv")
# * Sequences
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
