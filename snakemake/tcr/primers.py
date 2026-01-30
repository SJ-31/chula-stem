#!/usr/bin/env ipython

import json
from functools import partial
from pathlib import Path
from typing import Callable, Literal

import anndata as ad
import mudata as md
import polars as pl
import polars.selectors as cs
import pydna.tm as tm
import scirpy as ir
from pydna.design import Amplicon, pcr, primer_design
from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer

try:
    from snakemake.script import snakemake as smk
except ImportError:
    smk = type("snakemake", (), {"rule": None, "config": {None: None}, "log": [0]})

RCONFIG: dict = smk.config[smk.rule]

CHAINS: tuple = ("VJ_1", "VDJ_1")
WANTED_COLS = [
    "cdr3",
    "j_sequence_start",
    "j_sequence_end",
    "v_sequence_start",
    "v_call",
    "j_call",
    "d_call",
    "v_sequence_end",
    "sequence",
]


def get_cdr3_indices(
    airr: ad.AnnData,
    chain: str = "VJ_1",
    verify: bool = True,
) -> pl.DataFrame:
    pre_cdr3 = [
        "fwr1",
        "cdr1",
        "fwr2",
        "cdr2",
        "fwr3",
    ]
    tmp = pl.from_pandas(
        ir.get.airr(
            airr,
            ["sequence"] + pre_cdr3 + ["cdr3", "v_sequence_start"],
            chain=chain,
        ),
        include_index=True,
    ).rename(lambda x: x.removeprefix(f"{chain}_"))
    length_expr = [pl.col(c).str.len_chars() for c in pre_cdr3]
    start = f"{chain}_cdr3_start"
    df = (
        tmp.with_columns(*length_expr)
        .with_columns(
            (pl.col("v_sequence_start") + pl.sum_horizontal(pre_cdr3)).alias(start),
            pl.col("cdr3").str.len_chars().alias("len"),
        )
        .with_columns((pl.col(start) + pl.col("len")).alias(f"{chain}_cdr3_end"))
    )
    if verify:
        check = f"{chain}_cdr3_index_consistent"
        df = df.with_columns(
            pl.col("sequence")
            .str.slice(offset=pl.col(start), length=pl.col("len"))
            .alias(check)
        ).with_columns(pl.col(check) == pl.col("cdr3"))
    return df.select("index", cs.starts_with(chain))


def calc_tm(
    seq: str,
    method: Literal["neb", "default"] = "default",
    prodcode: str = "hstaq-0",
    **kws,
):
    """Wrapper function around pydna tm calculations
    See https://tmapi.neb.com/docs/productcodes for product codes
    """
    if method == "default":
        return tm.tm_default(seq, **kws)
    elif method == "neb":
        return tm.tm_neb(seq, prodcode=prodcode, **kws)


# TODO: add in support for adding arbitrary overhang sequences
def get_cdr3_fusion_primers(
    airr_data: dict,
    chain: str,
    how: Literal["pydna", "manual"] = "pydna",
    forward_len: int = 30,
    reverse_len: int = 30,
    tm_func: Callable[[str], float] = tm.tm_default,
    target_tm: float = 55.0,
) -> dict:
    full: Dseqrecord = Dseqrecord(airr_data[f"{chain}_sequence"])
    p1_start, p1_end = (
        airr_data[f"{chain}_v_sequence_start"],
        airr_data[f"{chain}_cdr3_end"],
    )
    p2_start, p2_end = (
        airr_data[f"{chain}_cdr3_start"],
        airr_data[f"{chain}_j_sequence_end"],
    )
    results = {}
    for i, (start, end) in enumerate([(p1_start, p1_end), (p2_start, p2_end)]):
        if how == "pydna":
            target = Dseqrecord(airr_data[f"{chain}_sequence"][start:end])
            amp = primer_design(target, target_tm=target_tm, tm_func=tm_func)
        elif how == "manual":
            fwd = full[start : start + forward_len]
            rev = full[end - reverse_len : end].reverse_complement()
            primers = (Primer(fwd, id=f"f{end-start}"), Primer(rev, id=f"r{end-start}"))
            amp = pcr(primers, full)
        results[f"pair {i+1}"] = amp

    return results


# * Formatting functions


class PrimerFmt:
    def __init__(self, how: Literal["text", "df", "json"] = "df", **kws) -> None:
        self.kws: dict = kws
        self.how: Literal["text", "df", "json"] = how

    def __call__(self, primers: dict[str, Amplicon]) -> str | dict | pl.DataFrame:
        if self.how == "text":
            return self._to_text(primers)
        elif self.how == "df":
            return self._to_df(primers)
        return self._2dict(primers)

    def _2dict(self, primers: dict[str, Amplicon]) -> dict:
        result = {}
        for name, amp in primers.items():
            result[name] = {"figure": amp.figure()}
            for p, pname in zip(amp.primers(), ("Forward", "Reverse")):
                if p is None:
                    continue
                seq = str(p.seq)
                result[pname] = {
                    "seq": seq,
                    "tm": calc_tm(seq, **self.kws),
                    "gc_content": p.gc(),
                    "length": len(seq),
                }
        return result

    def _to_text(self, primers: dict[str, Amplicon]) -> str:
        lines = []
        for name, amp in primers.items():
            lines.extend([f"Primer set: {name}"])
            fig = amp.figure()
            lines.extend([fig, ""])
            for p, pname in zip(amp.primers(), ("Forward", "Reverse")):
                if p is not None:
                    lines.append(pname)
                    seq = str(p.seq)
                    tm = calc_tm(seq, **self.kws)
                    lines.append(f"\tSequence: {seq}")
                    lines.append(
                        f"\tLength: {len(seq)}\tTM: {tm}\tGC content: {p.gc()}"
                    )
            lines.append("\n")
        return "\n".join(lines)

    def _to_df(
        self,
        primers: dict[str, Amplicon],
    ) -> str | pl.DataFrame:
        tmp = {
            "primer_set": [],
            "direction": [],
            "sequence": [],
            "tm": [],
            "gc_content": [],
        }
        for name, amp in primers.items():
            for p, pname in zip(amp.primers(), ("Forward", "Reverse")):
                if p is not None:
                    tmp["primer_set"].append(name)
                    tmp["direction"].append(pname)
                    seq = str(p.seq)
                    tm = calc_tm(seq, **self.kws)
                    tmp["sequence"].append(seq)
                    tmp["tm"].append(tm)
                    tmp["gc_content"].append(p.gc())
        return pl.DataFrame(tmp)


def format_one_sample(
    airr_data: dict, cfg: dict, tm_func=tm.tm_default
) -> tuple[pl.DataFrame, dict]:
    primer_dfs: list = []
    dct_result = {}
    for chain in CHAINS:
        primers = get_cdr3_fusion_primers(
            airr_data,
            chain=chain,
            tm_func=tm_func,
            target_tm=cfg.get("target_tm", 55),
            how=cfg.get("how", "pydna"),
        )
        fmt: pl.DataFrame = PrimerFmt("df")(primers)
        exprs = []
        dct_result[chain] = PrimerFmt("json")(primers)
        for key in [
            f"{chain}_sequence",
            f"{chain}_v_call",
            f"{chain}_d_call",
            f"{chain}_j_call",
        ]:
            new_key = key.removeprefix(f"{chain}_")
            val = airr_data.get(key)
            exprs.append(pl.lit(val).alias(new_key))
            dct_result[chain][new_key] = val
        fmt = fmt.with_columns(
            pl.lit(chain).alias("chain"),
            pl.lit(airr_data["index"]).alias("index"),
        ).rename({"sequence": "template_sequence"})
        primer_dfs.append(fmt)
    return pl.concat(primer_dfs, how="diagonal_relaxed"), dct_result


def main(airr: ad.AnnData, out_tabular: Path, out_json: Path, cfg: dict):
    meta_cols = ["Sample_Name", "chain_pairing", "clone_id"]
    fields = airr.obsm["airr"].fields
    df = (
        pl.from_pandas(
            ir.get.airr(airr, WANTED_COLS, chain=CHAINS),
            include_index=True,
        )
        .join(
            pl.from_pandas(airr.obs.loc[:, meta_cols], include_index=True),
            on="index",
        )
        .filter(pl.col("chain_pairing") == "single pair")
    )
    if "cdr3_start" not in fields and "cdr3_end" not in fields:
        for chain in CHAINS:
            df = df.join(
                get_cdr3_indices(airr, chain=chain, verify=True), on="index"
            ).filter(pl.col(f"{chain}_cdr3_index_consistent"))
    df = (
        df.filter(
            pl.struct("Sample_Name", "clone_id").map_elements(
                lambda x: [x["clone_id"], x["Sample_Name"]] in cfg["wanted_clones"],
                return_dtype=pl.Boolean,
            )
        )
        .group_by("Sample_Name", "clone_id")
        .agg(pl.all().first())
    )  # WARNING: taking just the first sequence is kinda arbitrary, maybe do something
    # better
    to_json_raw = []
    primer_dfs = []
    if tm_kws := cfg.get("tm_function"):
        tm_func = partial(calc_tm, **tm_kws)
    else:
        tm_func = tm.tm_default
    for row in df.iter_rows(named=True):
        cur, dct = format_one_sample(row, cfg=cfg, tm_func=tm_func)
        for col in meta_cols:
            dct[col] = row[col]
        to_json_raw.append(dct)
        primer_dfs.append(cur)
    final: pl.DataFrame = pl.concat(primer_dfs, how="diagonal_relaxed").join(
        df.select(*meta_cols, "index"), on="index"
    )
    with open(out_json, "w") as f:
        json.dump(to_json_raw, f)
    final.write_csv(out_tabular)


# * Entry point


def create_primers():
    airr = md.read_h5mu(smk.input[0])["airr"]
    main(
        airr=airr,
        cfg=RCONFIG,
        out_tabular=smk.output["primers_tabular"],
        out_json=smk.output["primers_json"],
    )


if fn := globals().get(smk.rule):
    fn()
