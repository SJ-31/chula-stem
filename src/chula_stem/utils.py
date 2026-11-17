#!/usr/bin/env ipython

from __future__ import annotations

import importlib.resources as res
from collections.abc import Sequence
from functools import reduce
from pathlib import Path
from typing import Callable, ClassVar, Literal, TypeAlias

import click
import polars as pl
from attrs import define


def str_split_unique(
    df: pl.DataFrame,
    col: str,
    separator: str,
    collapse: str | None = None,
    processing_fn: Callable = str.strip,
) -> pl.DataFrame:
    collapse = collapse if collapse is not None else separator
    return df.with_columns(
        pl.col(col)
        .str.split(separator)
        .list.unique()
        .map_elements(
            lambda str_list: list(map(lambda y: processing_fn(y), str_list)),
            return_dtype=pl.List(pl.String),
        )
        .list.join(collapse)
    )


def empty_string2null(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(
        pl.when(pl.col(pl.String).str.len_chars() == 0)
        .then(None)
        .otherwise(pl.col(pl.String))
        .name.keep()
    )


def contain_join(
    x, y, x_start, y_start, y_end, on=None, x_end=None, suffix="_right"
) -> pl.DataFrame:
    """
    Join x and y on rows where columns `x_start` and `x_end` (optional) are contained
        by the range of `y_start` and `y_end`
    Does not accept true overlapping joins
    """
    predicates = []
    if on:
        predicates.extend([pl.col(o) == pl.col(f"{o}{suffix}") for o in on])
    predicates.extend(
        [pl.col(y_start) <= pl.col(x_start), pl.col(x_start) <= pl.col(y_end)]
    )
    if x_end:
        predicates.extend([pl.col(x_end) <= pl.col(y_end)])
    return x.join_where(y, *predicates)


def add_loc(
    df: pl.DataFrame,
    chr_col: str = "chromosome",
    start_col: str = "start",
    end_col: str = "end",
    loc_name: str = "Locus",
) -> pl.DataFrame:
    return (
        df.with_columns(__range=pl.concat_str([start_col, end_col], separator="-"))
        .with_columns(
            pl.concat_str([chr_col, "__range"], separator=":").alias(loc_name)
        )
        .drop("__range")
    )


def parse_multiqc_vep(input: str, output: str = "") -> pl.DataFrame:
    vep = pl.read_csv(input, separator="\t")
    vep_data: dict = {"sample": [], "category": [], "key": [], "value": []}
    cols = [
        "Variant classes",
        "Consequences (most severe)",
        "Consequences (all)",
        "Coding consequences",
        "Variants by chromosome",
        "Position in protein",
        "General statistics",
    ]
    for row in vep.iter_rows(named=True):
        sample: str = row["Sample"]
        for col in cols:
            parsed = ast.literal_eval(row[col])
            vep_data["sample"].extend([sample] * len(parsed))
            vep_data["category"].extend([col] * len(parsed))
            vep_data["key"].extend(parsed.keys())
            vep_data["value"].extend(parsed.values())
    df = pl.DataFrame(vep_data)
    if output:
        df.write_csv(output, separator="\t")
    return df


@click.command()
@click.option(
    "-d",
    "--dir",
    default=".",
    help="Directory containing the count files.",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.option("-p", "--pattern", default=None, help="Glob pattern to match files.")
@click.option("-i", "--id-col", default=0, help="Index of the ID column.", type=int)
@click.option(
    "-c", "--count-col", default=1, help="Index of the count column.", type=int
)
@click.option(
    "-o",
    "--outfile",
    default="counts.csv",
    help="Output file name.",
    type=click.Path(writable=True),
)
@click.option("-s", "--sep", default="\t", help="Column separator.")
@click.option(
    "-H", "--has-header", is_flag=True, help="Specify if input files have a header row."
)
def combine_counts(
    dir=".",
    pattern=None,
    id_col=0,
    count_col=1,
    outfile="counts.csv",
    sep="\t",
    has_header=False,
) -> None:
    if pattern:
        files = Path(dir).glob(pattern)
    else:
        files = Path(dir).iterdir()

    def read_fn(file, name) -> pl.DataFrame:
        df = pl.read_csv(
            file,
            separator=sep,
            has_header=has_header,
            columns=[id_col, count_col],
            new_columns=["feature_id", name],
        )
        return df

    dfs: list[pl.DataFrame] = [read_fn(f, f.stem) for f in files]
    combined: pl.DataFrame = reduce(
        lambda x, y: x.join(y, on="feature_id", how="full")
        .with_columns(feature_id=pl.coalesce(["feature_id", "feature_id_right"]))
        .drop("feature_id_right"),
        dfs,
    )
    combined.write_csv(outfile, null_value="NA")


def do_call(fn, pars: dict | None):
    if pars is None:
        return fn(**{})
    else:
        return fn(**pars)


def read_existing(
    filename: Path, expr: Callable, read_fn: Callable | None = None
) -> None:
    if filename.exists() and read_fn is not None:
        return read_fn(filename)
    elif filename.exists():
        return
    else:
        return expr(filename)


# * Ensembl annotations

ID_FIELDS: TypeAlias = Literal[
    "ensembl",
    "ensembl_transcript",
    "mane_select",
    "mane_select_clinical",
    "gene_name",
    "hgnc_symbol",
    "ncbi",
    "ensembl_protein",
    "refseq_protein",
    "refseq_mrna",
    "refseq_ncrna",
]


@define
class Ensembl116:
    path: ClassVar[Path] = Path(res.files("chula_stem")) / "data" / "ensembl116"
    rename: ClassVar[dict] = {
        "Transcript stable ID": "ensembl_transcript",
        "Transcript stable ID version": "ensembl_transcript",
        "Chromosome/scaffold name": "chr",
        "RefSeq match transcript (MANE Select)": "mane_select",
        "RefSeq match transcript (MANE Plus Clinical)": "mane_select_clinical",
        "Gene stable ID": "ensembl",
        "Gene stable ID version": "ensembl",
        "Gene name": "gene_name",
        "Gene start (bp)": "start",
        "Gene end (bp)": "end",
        "Transcript start (bp)": "start",
        "Transcript end (bp)": "end",
        "NCBI gene (formerly Entrezgene) ID": "ncbi",
        "HGNC symbol": "hgnc_symbol",
    }

    @staticmethod
    def coords(
        what: Literal["transcript", "gene"] = "gene", versioned: bool = False
    ) -> pl.DataFrame:
        if what == "transcript":
            path = Ensembl116.path / "transcript_coords.parquet"
            col = "Transcript stable ID version"
        else:
            path = Ensembl116.path / "gene_coords.parquet"
            col = "Gene stable ID version"
        df = pl.read_parquet(path)
        if not versioned:
            df = df.drop(col).unique()
        else:
            df = df.drop(col.removesuffix(" version"))
        to_rename = {k: v for k, v in Ensembl116.rename.items() if k in df.columns}
        return df.rename(to_rename)

    @staticmethod
    def name_mapping(
        old: ID_FIELDS,
        new: ID_FIELDS = "ensembl_transcript",
        versioned: bool = False,
    ) -> dict:
        assert old != new
        remap = {
            "ensembl": "Gene stable ID version",
            "refseq_protein": "RefSeq peptide ID",
            "refseq_mrna": "RefSeq mRNA ID",
            "refseq_ncrna": "RefSeq ncRNA ID",
            "ensembl_protein": "Protein stable ID version",
            "ensembl_transcript": "Transcript stable ID version",
            "mane_select": "RefSeq match transcript (MANE Select)",
            "mane_select_clinical": "RefSeq match transcript (MANE Plus Clinical)",
            "hgnc_symbol": "HGNC symbol",
            "ncbi": "NCBI gene (formerly Entrezgene) ID",
            "gene_name": "Gene name",
        }
        oc, nc = remap[old], remap[new]
        refseq_id_file = Ensembl116.path / "ids_refseq.parquet"
        id_file = Ensembl116.path / "ids.parquet"
        join_refseq = "Transcript stable ID version"
        o_from_refseq: bool = old.startswith("refseq")
        n_from_refseq: bool = new.startswith("refseq")

        if (
            (o_from_refseq and n_from_refseq)
            or (n_from_refseq and oc == join_refseq)
            or (o_from_refseq and nc == join_refseq)
        ):
            lf = pl.scan_parquet(refseq_id_file)

        elif o_from_refseq or n_from_refseq:
            from_refseq = oc if o_from_refseq else nc
            from_main = nc if o_from_refseq else oc
            lf = (
                pl.scan_parquet(id_file)
                .select([from_main, join_refseq])
                .unique()
                .join(
                    pl.scan_parquet(refseq_id_file)
                    .select([from_refseq, join_refseq])
                    .unique(),
                    on=join_refseq,
                )
            )
        else:
            lf = pl.scan_parquet(id_file)
        df = (
            lf.select([oc, nc])
            .with_columns(
                pl.when(not versioned and oc.endswith("version"))
                .then(pl.col(oc).str.replace("\\..*$", ""))
                .otherwise(pl.col(oc)),
                pl.when(not versioned and nc.endswith("version"))
                .then(pl.col(nc).str.replace("\\..*$", ""))
                .otherwise(pl.col(nc)),
            )
            .unique()
            .filter((pl.col(oc).is_not_null()) & (pl.col(nc).is_not_null()))
            .collect()
        )
        if df[oc].is_duplicated().any():
            print(f"WARNING: multiple ids in `{new}` map to the same id in `{old}`")
        return dict(zip(df[oc], df[nc]))

    @staticmethod
    def remap_hgvs(
        df: pl.DataFrame,
        col: str,
        old: ID_FIELDS,
        new: ID_FIELDS,
        versioned: bool = False,
        keep_failed: bool = True,
    ) -> pl.DataFrame:
        mapping = Ensembl116.name_mapping(old, new, versioned=versioned)
        tmp = (
            df.with_row_index()
            .select([col, "index"])
            .with_columns(pl.col(col).str.split(":").list.to_struct(["id", "change"]))
            .unnest(col)
            .with_columns(
                pl.when(versioned)
                .then(pl.col("id"))
                .otherwise(pl.col("id").str.replace("\\..*$", ""))
            )
            .with_columns(
                pl.when(keep_failed)
                .then(pl.col("id").replace(mapping))
                .otherwise(pl.col("id").replace_strict(mapping, default=""))
                .alias("id")
            )
            .with_columns((pl.col("id") + ":" + pl.col("change")).alias(col))
            .drop(["id", "change"])
        )
        df = df.with_row_index().drop(col).join(tmp, on="index")
        if not keep_failed:
            before = df.height
            df = df.filter(~pl.col(col).str.starts_with(":"))
            after = df.height
            print(f"Remapping failed for {before-after} rows")
        return df

    @staticmethod
    def biomart_attrs() -> set:
        return set(
            (Ensembl116.path / "biomart_attributes.txt").read_text().splitlines()
        )


def construct_biomart_query(
    attributes: Sequence, dataset_name="hsapiens_gene_ensembl"
) -> tuple[list[str], str]:
    """
    Constructs a query for retrieving attributes from ensembl biomart

    """
    valid_all: set = Ensembl116.biomart_attrs()
    valid_attrs: list = []
    for attr in attributes:
        if attr not in valid_all:
            print(
                f"WARNING: requested attribute {attr} is not a valid biomart attribute"
            )
        else:
            valid_attrs.append(attr)
    query = f"""<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6"><Dataset name = "{dataset_name}" interface="default">
    """
    for va in valid_attrs:
        query = f'{query}<Attribute name = "{va}" />'
    query = f"{query}</Dataset></Query>"
    return valid_attrs, query
