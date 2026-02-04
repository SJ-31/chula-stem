#!/usr/bin/env ipython

import json
import re
import subprocess as sp
from collections.abc import Sequence
from functools import partial
from pathlib import Path
from typing import Callable, Literal

import anndata as ad
import mudata as md
import polars as pl
import polars.selectors as cs
import pydna.tm as tm
import pytest
import scirpy as ir
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord
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
    "d_sequence_start",
    "d_sequence_end",
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

    @property
    def chain_name(self):
        name = None
        for gene in ("c", "v", "d", "j"):
            if call := self.get(f"{gene}_call"):
                return get_chain_name_from_call(call)
        return name

    def __contains__(self, key):
        return f"{self.chain}_{key}" in self.data

    def start_end(self, key):
        return self.get(f"{key}_sequence_start"), self.get(f"{key}_sequence_end")

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
    if call.startswith("TRB"):
        return "TRB"
    if call.startswith("TRG"):
        return "TRG"
    if call.startswith("TRD") or re.match(".*TRA.*DV.*", call):
        return "TRD"
    if call.startswith("TRA"):
        return "TRA"
    raise ValueError(f"Chain name couldn't be determined from call {call}")


def get_imgt_seqs(
    file: Path,
    region: Literal["c", "d", "j", "v", "leader"],
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


def get_best_alignment(
    query: DNA, references: list[DNA], free_ends
) -> tuple[float, DNA]:
    cur_best = (
        pair_align(query, references[0], free_ends=free_ends).score,
        references[0],
    )
    for s in references[1:]:
        alignment = pair_align(query, s, free_ends=free_ends)
        if alignment.score > cur_best[0]:
            cur_best = (alignment.score, s)
    return cur_best


def infer_air_endpoint(
    cdata: ChainData,
    endpoint: Literal["c", "leader"] = "c",
    c_try_match_allele: bool = True,
) -> DNA:
    """
    Infer sequence of an end segment (leader or C gene) of a TCR from available data

    Parameters
    ----------
    c_try_match_allele : boolean
        When inferring the C gene, if the c call is generic and does not match an allele,
        then return the first allele unless this is false.
        Otherwise, the sequence of the allele with the highest alignment score is returned
    """
    data_dir: Path = Path(
        sp.run("stitchr -dd", shell=True, capture_output=True).stdout.decode().strip()
    )
    chain_name: str = cdata.chain_name
    seq_file = data_dir / "HUMAN" / f"{chain_name}.fasta"
    if not seq_file.exists():
        raise ValueError(f"{chain_name}.fasta not found in stitchr data directory")
    candidates = get_imgt_seqs(seq_file, region=endpoint)
    if not candidates:
        raise ValueError(
            f"No sequences for {endpoint} could be found in {seq_file}. Try checking its contents"
        )
    if (endpoint == "c" and not c_try_match_allele) or len(candidates) == 1:
        return candidates[0]
    try:
        full = cdata["sequence"]
    except KeyError:
        raise ValueError("Cannot infer endpoint without the full TCR sequence")
    if endpoint == "c":
        seq = full[cdata["j_sequence_end"] :]
        free_ends = [True, False]
    else:
        seq = full[: cdata["v_sequence_start"]]
        free_ends = [False, True]
    score, best_seq = get_best_alignment(DNA(seq), candidates, free_ends=free_ends)
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
        new_c_seq: DNA = infer_air_endpoint(cdata, "c", c_try_match_allele=True)
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


def filter_wanted_clones(airr: ad.AnnData, cfg: dict, extras=()):
    fields = airr.obsm["airr"].fields
    df = (
        pl.from_pandas(
            ir.get.airr(airr, AIRR_WANTED_COLS + list(extras), chain=CHAINS),
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


def get_seq_from_cfg(
    cfg: dict,
    key: tuple[str, str] | str,
    name: str = "",
    prefix: str = "",
    default: str = "",
) -> DNA:
    if isinstance(key, str):
        val = cfg.get(key)
    else:
        val = cfg.get(key[0], {}).get(key[1])
    val = val or default
    if (path := Path(val)).exists() and val != "":
        seq = DNA.read(path, "fasta")
        seq.metadata["id"] = f"{prefix}{seq.metadata['id']}"
        return seq
    return DNA(val, metadata={"id": name})


def get_construct_chain_regions(
    cdata: ChainData, cfg: dict, acc, gene_offset, allow_stop: bool
) -> tuple[list, list, int]:
    """Combine sequences from ChainData object `cdata` into a construct
    Intended to be called by `create_one_construct`

    Parameters
    ----------
    param : argument

    Returns
    -------
    tuple of list, list, int
    - The first list contains the ordered sequences in the construct
    [v_leader] <fwr1> <cdr1> <fwr2> <cdr2> <fwr3> <cdr3> <fwr4> [c_gene]
    - The second are construct features for visualization purposes, which includes
        all elements of the first list as well as the V, (D), J genes
    - Adjusted offset for incrementing intervals correctly in visualization
    """
    construct_seqs = []
    for_viz = []
    gene_offset -= cdata.get("v_sequence_start", 0)
    chain_name = get_chain_name_from_call(cdata["v_call"])
    if cfg.get("include_leader", False):
        leader = get_leader(cdata, chain_name, cfg)
        leader.metadata["interval"] = (acc, acc + len(leader))
        construct_seqs.append(leader)
        for_viz.append(leader)
        acc += len(leader)
        gene_offset += len(leader)
    for region in ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"):
        region_seq = DNA(cdata[region], metadata={"id": region})
        region_seq.metadata["interval"] = (acc, acc + len(region_seq))
        acc += len(region_seq)
        construct_seqs.append(region_seq)
        for_viz.append(region_seq)
    full = cdata.get("sequence")
    genes = ["v", "j"]
    if cdata.chain.startswith("VDJ"):
        genes.insert(1, "d")
    if full:
        add_genes_to_viz(full, for_viz, genes, cdata, gene_offset)
    if cfg.get("include_c", True):
        c_gene = get_c_gene(cdata, chain_name, cfg, allow_stop=allow_stop)
        c_gene.metadata["interval"] = (acc, acc + len(c_gene))
        acc += len(c_gene)
        construct_seqs.append(c_gene)
        for_viz.append(c_gene)
    return construct_seqs, for_viz, acc


def add_genes_to_viz(full, viz_list: list, genes, cdata, gene_offset):
    for g in genes:
        if not cdata.get(f"{g}_call"):
            continue
        start, end = cdata.start_end(g)
        seq = full[start:end]
        call = cdata[f"{g}_call"]
        meta = {
            "id": f"{g.upper()} gene: {call}",
            "interval": (start + gene_offset, end + gene_offset),
        }
        viz_list.append(DNA(seq, metadata=meta))


def get_leader(cdata, chain_name, cfg) -> DNA:
    leader = get_seq_from_cfg(
        cfg["sequences"], ("leader", chain_name), "V leader", "V leader"
    )
    if not leader:
        leader = infer_air_endpoint(cdata, "leader")
        leader_allele = leader.metadata["description"].split("|")[1]
        leader.metadata["id"] = f"V leader: {leader_allele}"
        logger.info(f"Inferred leader sequence as {leader_allele}")
    return leader


def get_c_gene(cdata, chain_name, cfg, allow_stop: bool = False) -> DNA:
    c_gene: DNA = get_seq_from_cfg(
        cfg.get("sequences", {}),
        ("c_gene", chain_name),
        "C gene",
        prefix="C gene: ",
        default="",
    )
    if not c_gene:
        c_gene = infer_air_endpoint(cdata, "c")
        c_allele = c_gene.metadata["description"].split("|")[1]
        c_gene.metadata["id"] = f"C gene: {c_allele}"
        logger.info(f"Inferred C gene as {c_allele}")
    last_codon = str(c_gene)[-3:]
    if last_codon in {"TAA", "TAG", "TGA"} and not allow_stop:
        c_gene = c_gene[:-3]
    return c_gene


def construct_add_flanking(
    cfg: dict,
    flank: Literal["five_prime", "three_prime"],
    construct_seqs: list,
    for_viz: list,
    acc: int,
) -> int:
    name = "5' flank" if flank == "five_prime" else "3' flank"
    flank_seq = get_seq_from_cfg(
        cfg.get("sequences", {}), ("flanking", flank), name, default=""
    )
    if flank_seq:
        flank_seq.metadata["interval"] = (acc, acc + len(flank_seq))
        construct_seqs.append(flank_seq)
        for_viz.append(flank_seq)
        acc += len(flank_seq)
    return acc


def plot_construct(
    sequence: str,
    features: list[DNA],
    cfg: dict,
    ignored: Sequence = ("cdr1", "cdr2", "fwr1", "fwr2", "fwr3", "fwr4"),
):
    graphic_features = []
    colormap = {
        "C gene": "#aaff32",
        "J gene": "#c04e01",
        "V gene": "#0485d1",
        "D gene": "#ffd1df",
        "V leader": "#bc13fe",
        "cdr3": "#ff028d",
        "linker": "#ff474c",
    }
    colormap.update(cfg.get("colormap", {}))
    for feature in features:
        start, end = feature.metadata["interval"]
        label = feature.metadata["id"]
        prefix = re.sub(":.*", "", label)
        color = colormap.get(prefix, "#fac205")
        if label not in ignored:
            gf = GraphicFeature(
                start=start + 1, end=end, label=label, strand=+1, color=color
            )
            graphic_features.append(gf)
    record = GraphicRecord(sequence=sequence, features=graphic_features)
    return record


def create_one_construct(airr_data: dict, chains, cfg: dict) -> DNA:
    tmp = []
    for_viz = []
    acc = 0
    seq_cfg: dict = cfg.get("sequences", {})
    gene_offset = acc = construct_add_flanking(cfg, "five_prime", tmp, for_viz, acc)
    for i, chain in enumerate(chains):
        is_final_chain = i == len(chains) - 1
        cdata = ChainData(chain, airr_data)
        cur_seqs, cur_viz, acc = get_construct_chain_regions(
            cdata, cfg, acc, gene_offset, allow_stop=is_final_chain
        )
        tmp.extend(cur_seqs)
        for_viz.extend(cur_viz)
        linker = get_seq_from_cfg(
            seq_cfg, "linker", "linker", prefix="linker: ", default=""
        )
        if linker and not is_final_chain:
            linker.metadata["interval"] = (acc, acc + len(linker))
            tmp.append(linker)
            for_viz.append(linker)
            acc += len(linker)
        gene_offset = acc
    construct_add_flanking(cfg, "three_prime", tmp, for_viz, acc)
    full_construct = "".join([str(seq) for seq in tmp])
    # for v in for_viz:
    #     logger.debug(f"{v.metadata['id']}, {v.metadata['interval']}")
    return full_construct, plot_construct(full_construct, for_viz, cfg)


def create_construct_main(
    airr: ad.AnnData, out_tabular: Path, out_json: Path, cfg: dict
):
    df = filter_wanted_clones(
        airr, cfg, extras=["fwr1", "fwr2", "fwr3", "fwr4", "cdr2", "cdr1"]
    )


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
