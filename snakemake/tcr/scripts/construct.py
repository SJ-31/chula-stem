#!/usr/bin/env python3
import re
import subprocess as sp
from collections.abc import Sequence
from dataclasses import asdict, dataclass, field
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


@beartype
def get_best_alignment(
    query: DNA, references: list[DNA], free_ends, return_alignment: bool = False
) -> tuple[float, DNA] | tuple[PairAlignResult, DNA]:
    best_seq = references[0]
    best_align = pair_align(query, best_seq, free_ends=free_ends)
    cur_best = (
        best_align.score,
        references[0],
    )
    for s in references[1:]:
        alignment = pair_align(query, s, free_ends=free_ends)
        if alignment.score > cur_best[0]:
            cur_best = (alignment.score, s)
            best_align = alignment
            best_seq = s
    if not return_alignment:
        return cur_best
    return best_align, best_seq


def plot_construct(
    sequence: str,
    features: list[DNA],
    cfg: dict,
    ignored: Sequence = ("CDR1", "CDR2", "FWR1", "FWR2", "FWR3", "FWR4"),
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
        seq: DNA = DNA.read(f, "fasta")
        seq.metadata["id"] = (
            f"{prefix}{seq.metadata['id']} {seq.metadata['description']}"
        )
        if metadata:
            seq.metadata.update(metadata)
        return seq

    sequences = []
    if isinstance(val, str):
        val = [val]
    for v in val:
        if (path := Path(v)).exists() and v:
            if path.is_file():
                sequences.append(read_one(path))
                continue
            elif path.is_dir():
                sequences.extend([read_one(p) for p in path.iterdir()])
        elif v:
            sequences.append(DNA(val, metadata={"id": name}))
    return sequences


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
                    metadata={"id": desc[1], "description": "|".join(desc)},
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

    Notes
    -----
    - For the C gene, this uses alignment between remaining sequence downstream of the J
    - For the V leader, uses alignment upstream of the V sequence
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
        free_ends = (True, False)
    else:
        seq = full[: cdata["v_sequence_start"]]
        free_ends = (False, True)
    query = DNA(seq)
    scored = [
        (pair_align(query, cand, free_ends=free_ends).score, cand)
        for cand in candidates
    ]
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
    check_has_start: bool | None = None
    # Either the leader sequence, or the start of the TCR block
    check_no_inframe_stop: bool = True
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
        self.ref_sequences: dict[tuple[str, str], dict[str, DNA]] | None = ref_sequences
        # Dictionary of (chain, gene) -> {allele->DNA}
        self.seen_stop_indices: set[int] = set()
        # Keep track of positions where stop codons were seen already
        self.ref_sequences_fastas: dict[str, Path] = {}
        if ref_seq_dirs := self.cfg.get("ref_seq_dirs"):
            for dir_or_file in (Path(d) for d in ref_seq_dirs):
                if dir_or_file.exists() and dir_or_file.is_dir():
                    for f in dir_or_file.iterdir():
                        self.ref_sequences_fastas[f.stem] = f
                elif dir_or_file.exists():
                    self.ref_sequences_fastas[dir_or_file.stem] = dir_or_file

    def _get_leaders(self, cdata, chain_name) -> list[DNA]:
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

    def _get_c_genes(self, cdata, chain_name, allow_stop: bool = False) -> list[DNA]:
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
                c.metadata["region_type"] = "c_gene"
                last_codon = str(c)[-3:]
                if last_codon in {"TAA", "TAG", "TGA"} and not allow_stop:
                    c = c[:-3]
                c_genes.append(c)
        return c_genes

    def _add_genes_to_viz(
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
            if g != "d":  # TODO: check why there are no d genes
                self._check_alignment(DNA(seq), g, call, chain=chain_name)
            meta = {
                "id": f"{g.upper()} gene: {call}",
                "interval": (start + gene_offset, end + gene_offset),
                "region_type": f"{g}_gene",
            }
            self.viz_data.append(DNA(seq, metadata=meta))

    def _add_chain_regions(
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
            leaders = self._get_leaders(cdata, chain_name)
            for l in leaders:
                l.metadata["interval"] = (acc, acc + len(l))
            leader = self._add_seq_check_inframe_stop(leaders, chain_index=chain_index)
            self.viz_data.append(leader)
            acc += len(leader)
            gene_offset += len(leader)
        for region in ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"):
            region_seq = DNA(
                cdata[region], metadata={"id": region.upper(), "region_type": region}
            )
            region_seq.metadata["interval"] = (acc, acc + len(region_seq))
            acc += len(region_seq)
            self._add_seq_check_inframe_stop([region_seq], chain_index)
            self.viz_data.append(region_seq)
        full = cdata.get("sequence")
        genes: list = ["v", "j"]
        if cdata.chain.startswith("VDJ"):
            genes.insert(1, "d")
        if full:
            self._add_genes_to_viz(
                full, genes, cdata, gene_offset, chain_name=chain_name
            )
        if self.cfg.get("include_c", True):
            c_genes = self._get_c_genes(cdata, chain_name, allow_stop=allow_stop)
            for c in c_genes:
                c.metadata["interval"] = (acc, acc + len(c))
            c_gene = self._add_seq_check_inframe_stop(
                c_genes, chain_index, allow_terminal=allow_stop
            )
            acc += len(c_gene)
            self.viz_data.append(c_gene)
        return acc

    def _add_flanking(
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
            added = self._add_seq_check_inframe_stop(seqs, 0)
            self.viz_data.append(added)
            acc += len(added)
        return acc

    @beartype
    def _check_alignment(
        self,
        sequence: DNA,
        gene: Literal["v", "j"],
        allele: str,
        chain: CHAIN_NAME,
    ):
        """
        Check the called V or J allele/gene with a reference sequence by alignment

        Notes
        -----
        IMGT reference sequences are obtained from stitchr's data directory, or the
            user-provided reference sequences in the config "ref_seq_dirs". If no
            match can be found from those sources, the best match will be returned
            from the IMGT sequences in stitchr for the given chain and gene

        TODO: create plots
        """
        if not self.check_with_reference:
            return
        if not self.ref_sequences:
            self.ref_sequences = {}
        if (chain, gene) not in self.ref_sequences:
            seq_file = get_stitchr_file(chain)
            seqs = get_imgt_seqs(seq_file, region=gene)
            self.ref_sequences[(chain, gene)] = {s.metadata["id"]: s for s in seqs}
        cur: dict[str, DNA] = self.ref_sequences[(chain, gene)]
        result = {}
        if match := cur.get(allele):
            aligned = pair_align(sequence, match)
        elif match_with_ext := self.ref_sequences_fastas.get(allele):
            tmp = next(SeqIO.parse(match_with_ext, "fasta"))
            match = DNA(
                str(tmp.seq),
                lowercase=True,
                metadata={"description": tmp.description, "id": tmp.id},
            )
            self.ref_sequences[(chain, gene)][allele] = match
            aligned = pair_align(sequence, match)
        else:
            logger.warning(
                f"The sequence {allele} couldn't be found in stitchr's IMGT database. Will return best match"
            )
            aligned, match = get_best_alignment(
                sequence, list(cur.values()), free_ends=True, return_alignment=True
            )
        assert isinstance(aligned, PairAlignResult)
        seq_aa = sequence.translate()
        match_aa = match.translate()
        aa_aligned = pair_align(seq_aa, match_aa)
        result["call"] = allele
        result["sequence_in_db"] = allele in cur
        result["matching_seq"] = match.metadata
        result["dna"] = {
            "score": aligned.score,
            "query_length": len(sequence),
            "alignment": aligned.paths[0].to_aligned((sequence, match)),
        }
        result["aa"] = {
            "score": aa_aligned.score,
            "query_length": len(seq_aa),
            "alignment": aa_aligned.paths[0].to_aligned((seq_aa, match_aa)),
        }
        self.vreport.gene_validation[gene] = result

    # def _plot_alignment_validation(self, vresult: dict):

    def _add_seq_check_inframe_stop(
        self,
        seqs: list[DNA],
        chain_index: int,
        allow_terminal: bool = False,
    ) -> DNA:
        """
        Append construct sequence `seq` to the internal list, and check if adding it would
        introduce stop codons

        Notes
        -----
        For the purposes of checking whether inframe stops are present, the 5' and 3'
        sequences are ignored. Only the TCR block is checked, following instruction of [1]

        References
        ----------
        [1] Afeyan AB, Wu CJ, Oliveira G. Rapid parallel reconstruction and specificity screening of hundreds of T cell receptors. Nat Protoc. 2025 Mar;20(3):539-586. doi: 10.1038/s41596-024-01061-4. Epub 2024 Nov 8. PMID: 39516267; PMCID: PMC11896752.
        """
        tracker: dict[int, bool] = {}
        stop_indices: dict[int, list] = {}

        def has_ifs(seq: DNA, idx: int) -> bool:
            block = DNA.concat([self.tcr_block] + [seq])
            stops: np.ndarray = block.translate().stops()
            has_terminal = stops[-1]
            indices = [
                i.item() for i in np.where(stops)[0] if i not in self.seen_stop_indices
            ]
            self.seen_stop_indices |= set(indices)
            stop_indices[idx] = indices
            if allow_terminal and has_terminal:
                res = len(indices) > 1
            else:
                res = len(indices) > 0
            tracker[idx] = res
            return res

        chosen_idx = 0
        for i, seq in enumerate(seqs):
            if not has_ifs(seq, i) or i == len(seqs) - 1:
                chosen_idx = i
                break
            logger.warning(
                f"Adding sequence {seq.metadata['id']} introduces stop. Trying another..."
            )
        chosen: DNA = seqs[chosen_idx]
        has_inframe_stops = tracker[chosen_idx]
        self.to_assemble.append(chosen)
        if self.vreport.check_no_inframe_stop and has_inframe_stops:
            self.vreport.check_no_inframe_stop = False
        if has_inframe_stops:
            offset = np.sum(
                [
                    len(s)
                    for s in self.to_assemble
                    if s.metadata["region_type"] == "five_prime"
                ]
            )
            indices = stop_indices[chosen_idx]
            meta = chosen.metadata.copy()
            meta["chain_number"] = chain_index
            meta["stop_codon_starts"] = [int(i * 3 + offset) for i in indices]
            self.vreport.inframe_stop_origins.append(meta)
        return chosen

    @property
    def tcr_block(self) -> DNA:
        return DNA.concat(
            [
                s
                for s in self.to_assemble
                if s.metadata["region_type"] not in {"five_prime", "three_prime"}
            ]
        )

    def _check_start_stop(self):
        """Check whether the start of the TCR block contains start and stop codons"""
        tcr_block: DNA = self.tcr_block
        stops: np.ndarray = tcr_block.translate().stops()
        self.vreport.check_has_terminal_stop = stops[-1].item()
        self.vreport.check_has_start = str(tcr_block[:3]) == "ATG"

    def assemble(self) -> dict:
        "Assemble construct from init data. The only function you need to call"
        acc = 0
        seq_cfg: dict = self.cfg.get("sequences", {})
        gene_offset = acc = self._add_flanking("five_prime", acc)
        for i, chain in enumerate(self.chains):
            is_final_chain = i == len(self.chains) - 1
            cdata = ChainData(chain, self.airr_data)
            acc = self._add_chain_regions(
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
                linker = self._add_seq_check_inframe_stop(linkers, chain_index=i)
                self.viz_data.append(linker)
                acc += len(linker)
            gene_offset = acc
        _ = self._add_flanking("three_prime", acc)
        full_construct = "".join([str(seq) for seq in self.to_assemble])
        with_names = {}
        self._check_start_stop()
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
    """
    Obtain CDR3 indices from airr data. Requires that the start of the V sequence is
    present (which is the equivalent to the start of fwr1).
    From that
    """
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
    """Subset airr object to the wanted clones specified in the config
    and generate necessary columns
    """
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
    reference_db = {}
    for row in df.iter_rows(named=True):
        prefix = f"cid{row["clone_id"]}_{row["Sample_Name"]}"
        cons = TCRConstruct(airr_data=row, chains=["VJ_1", "VDJ_1"], cfg=rconfig)
        result: dict = cons.assemble()
        if cons.ref_sequences:
            reference_db.update(cons.ref_sequences)
        to_fasta = [f">{prefix}_full_construct\n{result['sequence']}"]
        to_fasta.extend([f">{k}\n{str(v)}" for k, v in result["named"].items()])
        with open(outdir / f"{prefix}.fasta", "w") as f:
            f.write("\n".join(to_fasta))
        with open(outdir / f"{prefix}.yaml", "w") as f:
            yaml.safe_dump(asdict(result["validation"]), f, sort_keys=False)
        ax, _ = result["plot"].plot(figure_width=rconfig.get("figure_width", 4))
        ax.figure.savefig(
            outdir / f"{prefix}.png",
            dpi=rconfig.get("dpi", 500),
            bbox_inches="tight",
        )


if __name__ == "__main__":
    assemble_constructs()
