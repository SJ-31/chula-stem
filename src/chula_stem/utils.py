#!/usr/bin/env ipython

import ast
import importlib.resources as res
from collections.abc import Sequence
from functools import reduce
from pathlib import Path
from typing import Callable

import click
import polars as pl


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


def get_valid_biomart_attributes() -> set:
    return set((res.files("data") / "biomart_attributes.txt").read_text().splitlines())


def construct_biomart_query(
    attributes: Sequence, dataset_name="hsapiens_gene_ensembl"
) -> tuple[list[str], str]:
    """
    Constructs a query for retrieving attributes from ensembl biomart

    """
    valid_all: set = get_valid_biomart_attributes()
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
