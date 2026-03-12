#!/usr/bin/env python3
import re
import subprocess as sp
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, TypeAlias

import anndata as ad
import numpy as np
import polars as pl
import polars.selectors as cs
import scirpy as ir
import yaml
from beartype import beartype
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord
from loguru import logger
from skbio.alignment import pair_align
from skbio.alignment._pair import PairAlignResult
from skbio.metadata import IntervalMetadata
from skbio.sequence import DNA

CHAIN_TYPES: TypeAlias = Literal["VJ_1", "VDJ_1", "VJ_2", "VDJ_2"]
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
CHAINS: tuple = ("VJ_1", "VDJ_1")
META_COLS = ["Sample_Name", "chain_pairing", "clone_id"]

REGION_TYPES: TypeAlias = Literal[
    "fwr1",
    "cdr1",
    "fwr2",
    "cdr2",
    "fwr3",
    "cdr3",
    "fwr4",
    "five_prime",
    "three_prime",
    "leader",
    "linker",
    "v_gene",
    "c_gene",
    "d_gene",
    "j_gene",
]

CHAIN_NAME: TypeAlias = Literal["TRB", "TRA", "TRD", "TRG"]


class ChainData:
    "Helper class for retrieving data following AIRR standards from a specific chain of a single clone"

    def __init__(self, chain: str, data: dict) -> None:
        self.data: dict = data
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


def get_best_alignment(
    query: DNA, references: list[DNA], free_ends, return_alignment: bool = False
) -> tuple[float, DNA] | tuple[PairAlignResult, DNA]:
    best_align = pair_align(query, references[0], free_ends=free_ends)
    cur_best = (
        best_align.score,
        references[0],
    )
    for s in references[1:]:
        alignment = pair_align(query, s, free_ends=free_ends)
        if alignment.score > cur_best[0]:
            cur_best = (alignment.score, s)
            best_align = alignment
    if not return_alignment:
        return cur_best
    return best_align, s


def plot_construct(
    sequence: str,
    features: list[DNA],
    cfg: dict,
    ignored: Sequence = ("cdr1", "cdr2", "fwr1", "fwr2", "fwr3", "fwr4"),
):
    graphic_features = []
    colormap: dict[REGION_TYPES, str] = {
        "c_gene": "#aaff32",
        "j_gene": "#c04e01",
        "v_gene": "#0485d1",
        "d_gene": "#ffd1df",
        "leader": "#bc13fe",
        "cdr3": "#ff028d",
        "linker": "#ff474c",
    }
    colormap.update(cfg.get("colormap", {}) or {})
    for feature in features:
        start, end = feature.metadata["interval"]
        label = feature.metadata["id"]
        rt: REGION_TYPES = feature.metadata["region_type"]
        color = colormap.get(rt, "#fac205")
        if label not in ignored:
            gf = GraphicFeature(
                start=start + 1, end=end, label=label, strand=+1, color=color
            )
            graphic_features.append(gf)
    record = GraphicRecord(sequence=sequence, features=graphic_features)
    return record


def get_seq_from_cfg(
    cfg: dict,
    key: tuple[str, str] | str,
    name: str = "",
    prefix: str = "",
    default: str = "",
    metadata: dict | None = None,
) -> list[DNA]:
    if isinstance(key, str):
        val = cfg.get(key)
    else:
        val = cfg.get(key[0], {}).get(key[1])
    val = val or default

    def read_one(f: Path) -> DNA:
        seq: DNA = DNA.read(path, "fasta")
        seq.metadata["id"] = (
            f"{prefix}{seq.metadata['id']} {seq.metadata['description']}"
        )
        if metadata:
            seq.metadata.update(metadata)
        return seq

    if (path := Path(val)).exists() and val != "":
        if path.is_file():
            return [read_one(p) for p in path.iterdir()]
    if val:
        return [DNA(val, metadata={"id": name})]
    return []


def get_chain_name_from_call(call: str) -> CHAIN_NAME:
    "Calls are assumed to follow IMGT nomenclature"
    # https://imgt.org/IMGTScientificChart/Nomenclature/IMGTallelepolymorphism.html
    # WARNING: some alleles seem shared between TRA and TRD for the V and leaders
    # e.g. TRAV14/DV4*01
    for chain in ("TRB", "TRG", "TRD", "TRA"):
        if call.startswith(chain):
            return chain
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
                    seq,  # "id" is the allele name
                    metadata={"id": desc[0], "description": "|".join(desc)},
                )
            )
    return seqs


def air_endpoint_candidates(
    cdata: ChainData,
    endpoint: Literal["c", "leader"] = "c",
    c_try_match_allele: bool = True,
) -> list[DNA]:
    """
    Infer sequence of an end segment (leader or C gene) of a TCR from available data

    Parameters
    ----------
    c_try_match_allele : boolean
        When inferring the C gene, if the c call is generic and does not match an allele,
        then return the first allele unless this is false.
        Otherwise, the sequence of the allele with the highest alignment score is returned
    """
    chain_name: CHAIN_NAME = cdata.chain_name
    seq_file = get_stitchr_file(chain_name)
    candidates = get_imgt_seqs(seq_file, region=endpoint)
    if not candidates:
        raise ValueError(
            f"No sequences for {endpoint} could be found in {seq_file}. Try checking its contents"
        )
    if (endpoint == "c" and not c_try_match_allele) or len(candidates) == 1:
        return candidates
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
    query = DNA(seq)
    scored = [(pair_align(query, cand).score, cand) for cand in candidates]
    scored = [s[1] for s in sorted(scored, key=lambda x: x[0], reverse=True)]
    return scored


def get_stitchr_file(chain_name: CHAIN_NAME) -> Path:
    data_dir: Path = Path(
        sp.run("stitchr -dd", shell=True, capture_output=True).stdout.decode().strip()
    )
    seq_file = data_dir / "HUMAN" / f"{chain_name}.fasta"
    if not seq_file.exists():
        raise ValueError(f"{chain_name}.fasta not found in stitchr data directory")
    return seq_file


@dataclass
class ConstructValidation:
    check_has_terminal_stop: bool | None = None
    check_has_start: bool | None = (
        None  # Either the leader sequence, or the start of the TCR block
    )
    check_no_inframe_stop: bool | None = None
    inframe_stops: list = field(default_factory=list)
    inframe_stop_origins: list = field(default_factory=list)
    gene_validation: dict = field(default_factory=dict)
    # Sequences that, when added, produce an inframe stop codon


class TCRConstruct:
    """

    viz_data : list[DNA]
    Construct features for visualization purposes, which includes
            all elements of the first list as well as the V, (D), J genes
    """

    def __init__(
        self,
        airr_data: dict,
        chains: list[CHAIN_TYPES],
        cfg: dict,
        ref_sequences: dict | None = None,
    ) -> None:
        self.cfg: dict = cfg
        self.chains: list[CHAIN_TYPES] = chains
        self.check_with_reference = cfg.get("check_with_reference", True)
        self.with_leader: bool = cfg.get("include_leader", False)
        self.airr_data: dict = airr_data
        self.viz_data: list[DNA] = []
        self.to_assemble: list[DNA] = []
        self.vreport: ConstructValidation = ConstructValidation()
        self.viz_offset: int = 0
        self.ref_sequences: dict | None = ref_sequences

    def get_leaders(self, cdata, chain_name) -> list[DNA]:
        leaders = get_seq_from_cfg(
            self.cfg["sequences"], ("leader", chain_name), "V leader", "V leader"
        )
        if not leaders:
            leaders = air_endpoint_candidates(cdata, "leader")
            for l in leaders:
                leader_allele = l.metadata["description"].split("|")[1]
                l.metadata["id"] = f"V leader: {leader_allele}"
        if leaders:
            for l in leaders:
                l.metadata["region_type"] = "v_leader"
        return leaders

    def get_c_genes(self, cdata, chain_name, allow_stop: bool = False) -> list[DNA]:
        tmp: list[DNA] = get_seq_from_cfg(
            self.cfg.get("sequences", {}),
            ("c_gene", chain_name),
            "C gene",
            prefix="C gene: ",
            default="",
        )
        if not tmp:
            tmp = air_endpoint_candidates(cdata, "c")
            for c in tmp:
                c_allele = c.metadata["description"].split("|")[1]
                c.metadata["id"] = f"C gene: {c_allele}"
        c_genes = []
        if tmp:
            for c in tmp:
                c.metadata["id"] = "C gene"
                c.metadata["region_type"] = "c_gene"
                last_codon = str(c)[-3:]
                if last_codon in {"TAA", "TAG", "TGA"} and not allow_stop:
                    c = c[:-3]
                c_genes.append(c)
        return c_genes

    def add_genes_to_viz(
        self,
        full,
        genes: list[Literal["v", "j", "d"]],
        cdata: ChainData,
        gene_offset,
        chain_name: CHAIN_NAME,
    ):
        for g in genes:
            if not cdata.get(f"{g}_call"):
                continue
            start, end = cdata.start_end(g)
            seq = full[start:end]
            call: str = cdata[f"{g}_call"]
            self._check_alignment(seq, g, call, chain=chain_name)
            meta = {
                "id": f"{g.upper()} gene: {call}",
                "interval": (start + gene_offset, end + gene_offset),
            }
            self.viz_data.append(DNA(seq, metadata=meta))

    def add_chain_regions(
        self, cdata: ChainData, acc, gene_offset, allow_stop: bool, chain_index: int
    ) -> int:
        """Combine sequences from ChainData object `cdata` into a construct
        Intended to be called by `assemble_one_construct`
        TODO: add configuration for custom V, J sequences

        Returns
        -------
        tuple of list, list, int
        - The first list contains the ordered sequences in the construct
        [v_leader] <fwr1> <cdr1> <fwr2> <cdr2> <fwr3> <cdr3> <fwr4> [c_gene]
        - Adjusted offset for incrementing intervals correctly in visualization
        """
        gene_offset -= cdata.get("v_sequence_start", 0)
        chain_name: CHAIN_NAME = get_chain_name_from_call(cdata["j_call"])
        if self.with_leader:
            leaders = self.get_leaders(cdata, chain_name)
            for l in leaders:
                l.metadata["interval"] = (acc, acc + len(l))
            leader = self.add_seq_check_inframe_stop(leaders, chain_index=chain_index)
            self.viz_data.append(leader)
            acc += len(leader)
            gene_offset += len(leader)
        for region in ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"):
            region_seq = DNA(
                cdata[region], metadata={"id": region.upper(), "region_type": region}
            )
            region_seq.metadata["interval"] = (acc, acc + len(region_seq))
            acc += len(region_seq)
            self.add_seq_check_inframe_stop([region_seq], chain_index)
            self.viz_data.append(region_seq)
        full = cdata.get("sequence")
        genes: list = ["v", "j"]
        if cdata.chain.startswith("VDJ"):
            genes.insert(1, "d")
        if full:
            self.add_genes_to_viz(
                full, genes, cdata, gene_offset, chain_name=chain_name
            )
        if self.cfg.get("include_c", True):
            c_genes = self.get_c_genes(cdata, chain_name, allow_stop=allow_stop)
            for c in c_genes:
                c.metadata["interval"] = (acc, acc + len(c))
            c_gene = self.add_seq_check_inframe_stop(
                c_genes, chain_index, allow_terminal=allow_stop
            )
            acc += len(c_gene)
            self.viz_data.append(c_gene)
        return acc

    def add_flanking(
        self,
        flank: Literal["five_prime", "three_prime"],
        acc: int,
    ) -> int:
        name = "5' flank" if flank == "five_prime" else "3' flank"
        seqs = get_seq_from_cfg(
            self.cfg.get("sequences", {}), ("flanking", flank), name, default=""
        )
        if seqs:
            for flank_seq in seqs:
                flank_seq.metadata["interval"] = (acc, acc + len(flank_seq))
                flank_seq.metadata["region_type"] = flank
                flank_seq.metadata["id"] = name
            added = self.add_seq_check_inframe_stop(seqs, 0)
            self.viz_data.append(added)
            acc += len(added)
        return acc

    @beartype
    def _check_alignment(
        self,
        sequence: DNA,
        gene: Literal["v", "d", "j"],
        allele: str,
        chain: CHAIN_NAME,
    ):
        if not self.check_with_reference:
            return
        if not self.ref_sequences:
            self.ref_sequences = {}
        if gene not in self.ref_sequences:
            seq_file = get_stitchr_file(chain)
            seqs = get_imgt_seqs(seq_file, region=gene)
            self.ref_sequences[gene] = {s.metadata["id"]: s for s in seqs}
        cur: dict[str, DNA] = self.ref_sequences[gene]
        result = {}
        if match := cur.get(allele):
            aligned = pair_align(sequence, match)
        else:
            aligned, match = get_best_alignment(
                sequence, cur.values(), free_ends=True, return_alignment=True
            )
        assert isinstance(aligned, PairAlignResult)
        result["score"] = aligned.score
        result["query_length"] = len(sequence)
        result["matching_seq"] = match.metadata
        result["alignment"] = aligned.paths[0].to_aligned(sequence, match)
        self.vreport.gene_validation[gene] = result

    def add_seq_check_inframe_stop(
        self, seqs: list[DNA], chain_index: int, allow_terminal: bool = False
    ) -> DNA:
        """
        Append construct sequence `seq` to the internal list, and check if adding it would
        introduce stop codons
        """
        tracker: dict[int, bool] = {}
        stop_indices: dict[int, np.ndarray] = {}

        def has_ifs(seq: DNA, idx: int) -> bool:
            block = DNA.concat(self.to_assemble + [seq])
            stops: np.ndarray = block.translate().stops()
            n_stops = stops.sum().item()
            has_terminal = stops[-1] and allow_terminal
            res = n_stops > 1 if (allow_terminal and has_terminal) else n_stops > 0
            tracker[idx] = res
            stop_indices[idx] = stops
            return res

        block: DNA = DNA.concat(self.to_assemble)
        chosen_idx = 0
        for i, seq in enumerate(seqs):
            if not has_ifs(seq, i) or i == len(seqs) - 1:
                chosen_idx = i
                break
            logger.warning(
                f"Adding sequence {seq.metadata['id']} introduces stop. Trying another..."
            )
        has_inframe_stops = tracker[chosen_idx]
        self.to_assemble.append(seqs[chosen_idx])
        self.vreport.check_no_inframe_stop = not has_inframe_stops
        if has_inframe_stops:
            seq.metadata["chain_index"] = chain_index
            self.vreport.inframe_stop_origins.append(seq.metadata)
            self.vreport.inframe_stop_origins.extend(stop_indices[chosen_idx].tolist())
        return seqs[chosen_idx]

    def _check_start_stop(self):
        """Check whether the start of the TCR block contains start and stop codons"""
        tcr_block: DNA = DNA.concat(self.to_assemble)
        stops: np.ndarray = tcr_block.translate().stops()
        self.vreport.check_has_terminal_stop = stops[-1].item()
        self.vreport.check_has_start = str(tcr_block[:3]) == "ATG"

    def assemble(self) -> dict:
        acc = 0
        seq_cfg: dict = self.cfg.get("sequences", {})
        gene_offset = acc = self.add_flanking("five_prime", acc)
        for i, chain in enumerate(self.chains):
            is_final_chain = i == len(self.chains) - 1
            cdata = ChainData(chain, self.airr_data)
            acc = self.add_chain_regions(
                cdata, acc, gene_offset, allow_stop=is_final_chain, chain_index=i
            )
            linkers = get_seq_from_cfg(
                seq_cfg, "linker", "linker", prefix="linker: ", default=""
            )
            if linkers and not is_final_chain:
                for l in linkers:
                    l.metadata["interval"] = (acc, acc + len(l))
                    l.metadata["region_type"] = "linker"
                    l.metadata["id"] = "linker"
                linker = self.add_seq_check_inframe_stop(linkers, chain_index=i)
                self.viz_data.append(linker)
                acc += len(linker)
            gene_offset = acc
        _ = self.add_flanking("three_prime", acc)
        self._check_start_stop()
        full_construct = "".join([str(seq) for seq in self.to_assemble])
        with_names = {}
        for i, seq in enumerate(self.to_assemble):
            name = seq.metadata.get("id")
            seq.interval_metadata = IntervalMetadata(len(seq))
            seq.interval_metadata.add(bounds=[(0, len(seq))], metadata={"name": name})
            with_names[f"{i}_{name}"] = seq
        return {
            "sequence": full_construct,
            "plot": plot_construct(full_construct, self.viz_data, self.cfg),
            "named": with_names,
            "validation": self.vreport,
        }


# * Rule


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


def assemble_constructs():
    import mudata as md
    from snakemake.script import snakemake as smk

    airr = md.read_h5mu(smk.input[0])["airr"]
    outdir: Path = smk.params["outdir"]
    rconfig = smk.config[smk.rule]
    outdir.mkdir(exist_ok=True)
    df = filter_wanted_clones(
        airr,
        rconfig,
        extras=["fwr1", "fwr2", "fwr3", "fwr4", "cdr2", "cdr1"],
    )
    # TODO: support providing multiple C gene candidates in the config
    for row in df.iter_rows(named=True):
        prefix = f"cid{row["clone_id"]}_{row["Sample_Name"]}"
        cons = TCRConstruct(airr_data=row, chains=["VJ_1", "VDJ_1"], cfg=rconfig)
        result: dict = cons.assemble()
        to_fasta = [f">{prefix}_full_construct\n{result['sequence']}"]
        to_fasta.extend([f">{k}\n{str(v)}" for k, v in result["named"].items()])
        with open(outdir / f"{prefix}.fasta", "w") as f:
            f.write("\n".join(to_fasta))
        with open(outdir / f"{prefix}.yaml", "w") as f:
            yaml.safe_dump(result["validation"], f)
        ax, _ = result["plot"].plot(figure_width=rconfig.get("figure_width", 4))
        ax.figure.savefig(
            outdir / f"{prefix}.png",
            dpi=rconfig.get("dpi", 500),
            bbox_inches="tight",
        )
