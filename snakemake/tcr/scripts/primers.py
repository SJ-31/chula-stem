#!/usr/bin/env ipython

import json
import re
import subprocess as sp
from functools import partial
from pathlib import Path
from typing import Callable, Literal

import anndata as ad
import mudata as md
import numpy as np
import polars as pl
import polars.selectors as cs
import pydna.tm as tm
import pytest
import scirpy as ir
from Bio import SeqIO
from loguru import logger
from pydna.design import Amplicon, pcr, primer_design
from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer
from pyhere import here
from skbio.alignment import pair_align
from skbio.sequence import DNA

try:
    from snakemake.script import snakemake as smk
except ImportError:
    smk = type("snakemake", (), {"rule": None, "config": {None: None}, "log": [0]})

RCONFIG: dict = smk.config[smk.rule]

CHAINS: tuple = ("VJ_1", "VDJ_1")
AIRR_WANTED_COLS = [
    "cdr3",
    "j_sequence_start",
    "j_sequence_end",
    "v_sequence_start",
    "v_call",
    "j_call",
    "c_call",
    "d_call",
    "v_sequence_end",
    "sequence",
]
META_COLS = ["Sample_Name", "chain_pairing", "clone_id"]


# * Utilities


class ChainData:
    "Helper class for retrieving data following AIRR standards from a specific chain of a single clone"

    def __init__(self, chain: str, data: dict) -> None:
        self.data: str = data
        self.chain: str = chain

    def __contains__(self, key):
        return f"{self.chain}_{key}" in self.data

    def get(self, key, default=None):
        return self.data.get(f"{self.chain}_{key}", default)

    def __getitem__(self, key: str):
        return self.data[f"{self.chain}_{key}"]

    def __setitem__(self, key, val) -> None:
        self.data[f"{self.chain}_{key}"] = val


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


def get_chain_name_from_call(call: str) -> Literal["TRB", "TRA", "TRD", "TRG"]:
    "Calls are assumed to follow IMGT nomenclature"
    # https://imgt.org/IMGTScientificChart/Nomenclature/IMGTallelepolymorphism.html
    if re.match("^TRB", call):
        return "TRB"
    if re.match("^TRG", call):
        return "TRG"
    if re.match("^TRD", call) or re.match(".*TRA.*DV.*", call):
        return "TRD"
    if re.match("^TRA", call):
        return "TRA"
    raise ValueError("Chain name couldn't be determined  from call")


def get_imgt_seqs(
    file: Path,
    region: Literal["c", "d", "j", "v", ""],
    species: str = "Homo sapiens",
) -> list[DNA]:
    """Return fasta entries corresponding to `gene` in a fasta file with IMGT headers"""
    seqs = []
    gene_key = f"{region.upper()}-REGION"
    for record in SeqIO.parse(file, "fasta"):
        desc = record.description.split("|")
        if (
            (desc[-1] == "~CONSTANT" and region == "c")
            or (desc[-1] == "~LEADER" and region == "leader")
            or (desc[4] == gene_key)
        ):
            seq = str(record.seq)
            if seq.islower():
                seq = seq.upper()
            seqs.append(
                DNA(
                    seq,
                    metadata={"id": desc[0], "description": "|".join(desc)},
                )
            )
    return seqs


def get_best_alignment(query: DNA, references: list[DNA]) -> tuple[float, DNA]:
    ends = (
        [True, False] if len(query) <= np.mean([len(s) for s in references]) else None
    )
    cur_best = (pair_align(query, references[0], free_ends=ends).score, references[0])
    for s in references[1:]:
        alignment = pair_align(query, s, free_ends=ends)
        if alignment.score > cur_best[0]:
            cur_best = (alignment.score, s)
    return cur_best


def get_c_sequence_from_call(cdata: ChainData, try_match_allele: bool = True) -> DNA:
    """
    Return the sequence of an AIR gene from its call (predicted identity)

    Parameters
    ----------
    try_match_allele : boolean
        If the call is generic and does not match an allele, then return the first allele
        for the gene sequence unless `try_match_allele` is True.
        If it is, the sequence of the allele with the highest alignment score is returned
    """
    call = cdata["c_call"]
    data_dir: Path = Path(
        sp.run("stitchr -dd", shell=True, capture_output=True).stdout.decode().strip()
    )
    chain_name: str = get_chain_name_from_call(call)
    seq_file = data_dir / "HUMAN" / f"{chain_name}.fasta"
    if not seq_file.exists():
        raise ValueError(f"{chain_name}.fasta not found in stitchr data directory")
    candidates = get_imgt_seqs(seq_file, "c")
    if not candidates:
        raise ValueError(
            f"No sequences for C gene could be found in {seq_file}. Try checking its contents"
        )
    if not try_match_allele or len(candidates) == 1:
        return candidates[0]
    full = cdata["sequence"]
    j_end = cdata["j_sequence_end"]
    seq = full[j_end + 1 :]
    score, best_seq = get_best_alignment(DNA(seq), candidates)
    return best_seq


def add_new_c(cdata: ChainData, c_seq: str) -> str:
    """
    Return full-length AIR sequence modified to have a new C gene sequence

    This function simply strips any sequence downstream of airr_data["j_sequence_end"]
    and replaces it with `c_seq`
    """
    j_end = cdata["j_sequence_end"]
    cut = cdata["sequence"][:j_end]
    return cut + c_seq


# TODO: add in support for adding arbitrary overhang sequences
def get_cdr3_fusion_primers(
    cdata: ChainData,
    how: Literal["pydna", "manual"] = "pydna",
    forward_len: int = 30,
    reverse_len: int = 30,
    end_gene: Literal["c", "j"] = "c",
    tm_func: Callable[[str], float] = tm.tm_default,
    target_tm: float = 55.0,
) -> tuple[dict, ChainData]:
    """
    Generate fusion primer pairs targeting the CDR3 region of a TCR in `airr_data`

    Parameters
    ----------
    end_gene : j, c
        Gene for which the reverse primer of the second primer pair should extend to
        Specifying j produces a fusion product of V-(D)-J; specifying c yields V-(D)-J-C

        In the case when data for the C gene is unavailable, a match is attempted using
        IMGT reference genes (in the stitchr data directory) and the (possibly partial)
        C gene sequence using the J sequence end

    Returns
    -------
    Dictionary of primer pairs and modified airr data dictionary

    """
    full: Dseqrecord = Dseqrecord(cdata["sequence"])
    end_gene = end_gene.lower()
    p1_start, p1_end = (cdata["v_sequence_start"], cdata["cdr3_end"])
    end_lookup: str = f"{end_gene}_sequence_end"
    full_seq: str = cdata["sequence"]
    if end_lookup not in cdata and end_gene == "c":
        logger.info("C gene sequence missing, retrieving from call...")
        new_c_seq: DNA = get_c_sequence_from_call(cdata, try_match_allele=True)
        full_seq = add_new_c(cdata, str(new_c_seq))
        cdata["c_call"] = new_c_seq.metadata["description"].split("|")[1]
        cdata["c_sequence_start"] = cdata["j_sequence_start"] + 1
        cdata["c_sequence_end"] = len(full_seq)
    elif end_lookup not in cdata:
        # TODO: you could also do this for the j call...
        raise ValueError("J gene sequence end not present in AIR data")
    p2_start, p2_end = (
        cdata["cdr3_start"],
        cdata[f"{end_gene}_sequence_end"],
    )
    results = {}
    for i, (start, end) in enumerate([(p1_start, p1_end), (p2_start, p2_end)]):
        if how == "pydna":
            target = Dseqrecord(full_seq[start:end])
            amp = primer_design(target, target_tm=target_tm, tm_func=tm_func)
        elif how == "manual":
            fwd = full[start : start + forward_len]
            rev = full[end - reverse_len : end].reverse_complement()
            primers = (Primer(fwd, id=f"f{end-start}"), Primer(rev, id=f"r{end-start}"))
            amp = pcr(primers, full)
        results[f"pair {i+1}"] = amp
    return results, cdata


# * Formatting functions


class PrimerFmt:
    def __init__(self, how: Literal["text", "df", "dict"] = "df", **kws) -> None:
        self.kws: dict = kws
        self.how: Literal["text", "df", "dict"] = how

    def __call__(self, primers: dict[str, Amplicon]) -> str | dict | pl.DataFrame:
        if self.how == "text":
            return self._to_text(primers)
        elif self.how == "df":
            return self._to_df(primers)
        return self._to_dict(primers)

    def _to_dict(self, primers: dict[str, Amplicon]) -> dict:
        result = {}
        for name, amp in primers.items():
            result[name] = {
                "figure": amp.figure(),
                "product_length": amp.accession.removesuffix("bp"),
            }
            for p, pname in zip(amp.primers(), ("Forward", "Reverse")):
                if p is None:
                    continue
                seq = str(p.seq)
                result[name][pname] = {
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
        cdata: ChainData = ChainData(chain, airr_data)
        primers, airr_update = get_cdr3_fusion_primers(
            cdata,
            tm_func=tm_func,
            target_tm=cfg.get("target_tm", 55),
            how=cfg.get("how", "pydna"),
            forward_len=cfg.get("forward_len", 30),
            reverse_len=cfg.get("reverse_len", 30),
            end_gene=cfg.get("end_gene", "c"),
        )
        fmt: pl.DataFrame = PrimerFmt("df")(primers)
        exprs = []
        dct_result[chain] = PrimerFmt("dict")(primers)
        for key in ["sequence", "v_call", "d_call", "j_call", "c_call"]:
            val = airr_update.get(key)
            exprs.append(pl.lit(val).alias(key))
            dct_result[chain][key] = val
        fmt = fmt.with_columns(
            pl.lit(chain).alias("chain"),
            pl.lit(airr_data["index"]).alias("index"),
        ).rename({"sequence": "template_sequence"})
        primer_dfs.append(fmt)
    return pl.concat(primer_dfs, how="diagonal_relaxed"), dct_result


def filter_wanted_clones(airr: ad.AnnData, cfg: dict):
    fields = airr.obsm["airr"].fields
    df = (
        pl.from_pandas(
            ir.get.airr(airr, AIRR_WANTED_COLS, chain=CHAINS),
            include_index=True,
        )
        .join(
            pl.from_pandas(airr.obs.loc[:, META_COLS], include_index=True),
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
    return df


# * Rules

# ** Construct creation


def create_one_construct(airr_data: dict, chain: str, cfg: str) -> DNA:
    v_seq = airr_data[f"{chain}"]


# def create_construct_main(airr: ad.AnnData, out_tabular: Path, out_json: Path, cfg: dict):

# ** Primer creation


def primers_main(airr: ad.AnnData, out_tabular: Path, out_json: Path, cfg: dict):
    df = filter_wanted_clones(airr, cfg)
    to_json_raw = []
    primer_dfs = []
    if tm_kws := cfg.get("tm_function"):
        tm_func = partial(calc_tm, **tm_kws)
    else:
        tm_func = tm.tm_default
    for row in df.iter_rows(named=True):
        cur, dct = format_one_sample(row, cfg=cfg, tm_func=tm_func)
        for col in META_COLS:
            dct[col] = row[col]
        to_json_raw.append(dct)
        primer_dfs.append(cur)
    final: pl.DataFrame = pl.concat(primer_dfs, how="diagonal_relaxed").join(
        df.select(*META_COLS, "index"), on="index"
    )
    with open(out_json, "w") as f:
        json.dump(to_json_raw, f)
    final.write_csv(out_tabular)


def create_primers():
    airr = md.read_h5mu(smk.input[0])["airr"]
    primers_main(
        airr=airr,
        cfg=RCONFIG,
        out_tabular=smk.output["primers_tabular"],
        out_json=smk.output["primers_json"],
    )


# * Entry point & tests

if fn := globals().get(smk.rule):
    fn()


@pytest.fixture
def airr_data_test():
    combined = md.read_h5mu(here("tests", "data", "airr_test.h5mu"))
    return combined["airr"]


def test_get_c(airr_data_test):
    chain = "VJ_1"
    df = pl.from_pandas(
        ir.get.airr(airr_data_test, AIRR_WANTED_COLS, chain=chain),
        include_index=True,
    )
    dct = next(df.iter_rows(named=True))
    full = dct[f"{chain}_sequence"]
    j_end = f"{chain}_j_sequence_end"
    j_start = f"{chain}_j_sequence_start"
    old_j = full[dct[j_start] : dct[j_end]]
    old_c = full[dct[j_end] :]
    new_c = get_c_sequence_from_call(dct, chain, try_match_allele=True)
    print(f"Old: {old_c}\nNew: {new_c}")
    print(f"Old full: {full}")
    new_full = add_new_c(dct, chain, new_c)
    print(f"New full: {new_full}")
    assert old_j == new_full[dct[j_start] : dct[j_end]]
