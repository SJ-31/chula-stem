#!/usr/bin/env ipython

from collections.abc import Sequence
from itertools import chain
from pathlib import Path
from subprocess import run
from tempfile import NamedTemporaryFile
from typing import Literal

import polars as pl
import rpy2.robjects as ro
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from matplotlib.figure import Figure
from pymsaviz import MsaViz

try:
    from snakemake.script import snakemake as smk
except ImportError:
    smk = type("snakemake", (), {"rule": None, "config": {}})


SEQ_COLS = ["sequence", "cdr1", "cdr2", "cdr3", "fwr1", "fwr2", "fwr3", "fwr4"]
VDJ_COLORS = {"v": "red", "j": "blue", "c": "green", "d": "orange", "d2": "purple"}


def get_stitchr_seqs(
    allele_names,
    species: str = "HUMAN",
    seqtype: Literal["~LEADER", "~CONSTANT", "~JOINING"] = "~LEADER",
) -> pl.DataFrame:
    """Return a dataframe containing the leader sequences of the specified alleles
    Requires stitchr to be installed
    """
    allele_set: set = set(allele_names)
    if len(allele_set) != len(allele_names):
        raise ValueError("WARNING: given list of allele names is not unique!")
    stitchr_dir = Path(
        run("stitchr -dd", shell=True, capture_output=True).stdout.decode().strip()
    )
    stitchr_dir = stitchr_dir / species
    fasta_to_load = []
    genes = {"TRA", "TRB", "TRD", "TRG"}
    for allele in allele_names:
        if not allele:
            continue
        seen_gene = ""
        for g in genes:
            if allele.startswith(g):
                gene_file = stitchr_dir / f"{g}.fasta"
                if not gene_file.exists():
                    raise ValueError(
                        f"gene file {g}.fasta is missing from the stitchr directory"
                    )
                seen_gene = g
                fasta_to_load.append(gene_file)
                break
        if not genes:
            break
        if seen_gene:
            genes.remove(seen_gene)

    seen_alleles = set()
    seq_dict = {"allele": [], "sequence": []}
    for seqrecord in chain(*(SeqIO.parse(fa, "fasta") for fa in fasta_to_load)):
        seqrecord: SeqRecord
        splits = seqrecord.description.split("|")
        if not splits or len(splits) < 2:
            continue
        allele, feature_type = splits[1], splits[-1]
        if (
            allele in allele_set
            and feature_type == seqtype
            and allele not in seen_alleles
        ):
            seen_alleles.add(allele)  # Precaution
            seq_dict["allele"].append(allele)
            seq_dict["sequence"].append(str(seqrecord.seq).upper())
        if len(seen_alleles) == len(allele_set):
            break
    return pl.DataFrame(seq_dict)


def global_local_alignment(
    pattern: dict[str, str],
    subject: str,
    outfile: Path | None = None,
    subject_name: str = "sequence",
    **kws,
) -> dict:
    """Wrapper for global-local alignment as implemented in pwalign"""
    ro.r("library(Biostrings)")
    ro.globalenv["pattern"] = ro.ListVector(pattern)
    ro.globalenv["subject"] = subject
    ro.globalenv["subject_name"] = subject_name
    ro.r("pattern <- unlist(pattern)")
    b_pattern = ro.r("BStringSet(pattern)")
    b_subject = ro.r("BStringSet(subject)")
    kws.update({"pattern": b_pattern, "subject": b_subject, "type": "global-local"})
    ro.globalenv["kws"] = ro.ListVector(kws)
    ro.r("align <- do.call(pwalign::pairwiseAlignment, kws)")
    ro.r("aligned_patterns <- pwalign::aligned(align)")
    ro.r("names(aligned_patterns) <- names(pattern)")
    ro.r("names(kws$subject) <- subject_name")
    ro.r("string_set <- c(kws$subject, aligned_patterns)")
    if outfile is not None:
        ro.globalenv["out"] = str(outfile)
        ro.r("writeXStringSet(string_set, out)")
    as_char = ro.r("as.character(string_set)")
    seqs = {as_char.names[i]: as_char[i] for i in range(len(as_char))}
    return seqs


def align_vdj(
    df: pl.DataFrame,
    chain: str,
    outdir: Path | str,
    nucleotide: bool = True,
    id_col: str = "sequence_id",
    additional_seqs: dict | None = None,
) -> pl.DataFrame:
    """Wrapper function for producing global-local alignments of TCR segments onto
        the full-length sequence

    Parameters
    ----------
    df : pl.DataFrame containing AIR data e.g. obtained with ir.get.airr
    additional_seqs : dict | None
        Optional dictionary describing columns in `df` to also include in the alignment
        format is sequence_name->column_name
    """
    additional_seqs = additional_seqs or {}

    def align_one(id, row: dict):
        outfile = outdir.joinpath(id).with_suffix(".fasta")
        sequence: str = ""
        patterns = {}
        for col in SEQ_COLS:
            key = f"{chain}_{col}" if nucleotide else f"{chain}_{col}_aa"
            val = row.get(key)
            if val and col != "sequence":
                patterns[col] = val
            elif val:
                sequence = val
        for pat, col in additional_seqs.items():
            val = row.get(col)
            if val:
                patterns[pat] = val
        if len(patterns) <= 1 or not sequence:
            return False
        else:
            f = NamedTemporaryFile("w+t", suffix=".fasta")
        _ = global_local_alignment(
            pattern=patterns, subject=sequence, subject_name="full", outfile=outfile
        )
        if f is not None:
            f.close()
        return False

    tracker = {id_col: [], "alignment_success": []}
    outdir = Path(outdir) if not isinstance(outdir, Path) else outdir
    outdir.mkdir(exist_ok=True)
    for id, row in zip(df[id_col], df.iter_rows(named=True)):
        if id is not None:
            success = align_one(id, row)
            tracker[id_col].append(id)
            tracker["alignment_success"].append(success)

    return pl.DataFrame(tracker)


def plot_vdj(
    id,
    chain: str,
    df: pl.DataFrame,
    file: str | Path,
    id_col: str = "sequence_id",
    title_spec: Sequence[tuple[str, str]] = (),
    title_kws: dict | None = None,
    **kws,
) -> Figure:
    """Plot vdj with pymsaviz

    Parameters
    ----------
    df : pl.DataFrame
        DataFrame containing AIR data, namely gene segment calls and sequence boundaries
    chain : str
        AIR chain to plot, e.g. VJ_1, VDJ_1
    title_spec : Sequence[tuple[str, str]]
        Specify data to include in the plot title
        Each element of the sequence is a tuple of (label, column_name)
        The title will be formatted in order of the given elements,
           with 'label: column_name_value', or just 'column_name_value' if label is
           the empty string
    title_kws : dict |  None
        Optional dict of keywords passed to axis.title

    Returns
    -------
    matplotlib figure generated by pymsaviz
    """
    title_kws = title_kws or {
        "loc": "left",
        "fontsize": "large",
        "pad": 20,
        "weight": "bold",
    }
    row = df.filter(pl.col(id_col) == id).row(0, named=True)
    mv: MsaViz = MsaViz(file, **kws)
    for seq in ["d", "d2", "v", "j", "c"]:
        start = row.get(f"{chain}_{seq}_sequence_start")
        call = row.get(f"{chain}_{seq}_call")
        end = row.get(f"{chain}_{seq}_sequence_end")
        if start is None:
            continue
        anno = f"{seq.upper()} call: {call}"
        rc = VDJ_COLORS.get(seq, "black")
        mv.add_text_annotation((start, end), text=anno, range_color=rc)
    fig: Figure = mv.plotfig()
    tmp = []
    for label, col in title_spec:
        lab = f"{label}: " if label else ""
        if val := row.get(col):
            tmp.append(f"{lab}{val}")
    fig.axes[0].set_title("  ".join(tmp), **title_kws)
    return fig
