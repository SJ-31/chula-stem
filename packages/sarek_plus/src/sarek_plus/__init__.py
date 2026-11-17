#!/usr/bin/env python3

import json
from collections.abc import Callable
from pathlib import Path

import click

from sarek_plus.databases import get_therapy_df
from sarek_plus.report.main import VariantCallingReport


@click.command()
@click.option("-o", "--output", required=False, help="Name of output pdf")
@click.option(
    "-t",
    "--report_type",
    required=True,
    default=None,
    help="Type of report to be generated",
)
@click.option(
    "-s", "--specification", required=True, help="Path specification file (YAML)"
)
@click.option(
    "-d",
    "--tmpdir",
    required=False,
    default="./report_tmp",
    help="Path to temporary directory for report files",
)
def entry_point(report_type: str, specification: str, output: str, tmpdir: str) -> None:
    with open(specification, "r") as f:
        spec: dict = json.load(f)
    try:
        paths: dict = spec["paths"]
        assert isinstance(paths, dict)
    except (KeyError, AssertionError):
        raise ValueError("`paths` map must be specified!")
    out = output or paths["output"]
    kws = spec.get("misc", {})
    report = None
    if report_type == "variant_calling":
        report = VariantCallingReport(
            filename=out,
            metadata=spec.get("meta", {}),
            other_text=spec.get("other_text", {}),
            tmpdir=tmpdir,
            paths=spec["paths"],
            **kws,
        )
    if report:
        report.build()
        report.merge()


@click.command()
@click.option(
    "-i",
    "--input-files",
    required=True,
    help="Path to directory containing input files",
)
@click.option(
    "-c",
    "--gene-col",
    default="Gene",
    show_default=True,
    help="For dataframe input, the column containing gene names",
)
@click.option(
    "-f",
    "--filter-confident",
    is_flag=True,
    default=False,
    show_default=True,
    help="Whether or not to filter confident therapy entries",
)
@click.option(
    "-n",
    "--not-in-db",
    default="not_in_databases.txt",
    show_default=True,
    help="File to write genes that weren't found",
)
@click.option(
    "-i",
    "--ignore_failed_file",
    required=True,
    help="If supplied, will ignore all genes in this file",
)
@click.option(
    "-e", "--cache", required=True, help="DB file to cache data", default=None
)
@click.option(
    "-d",
    "--databases",
    required=True,
    help="Comma-delimited string specifying which databases to include in the lookup",
    default=None,
    type=str,
)
def therapy_cache_from_files(
    input_files: str,
    not_in_db: str,
    filter_confident: bool,
    cache: str,
    gene_col: str = "Genes",
):
    import polars as pl

    genes: set = set()
    for path in Path(input_files).rglob("*"):
        path: Path
        if path.is_dir():
            continue
        if path.suffix == ".tsv":
            df = pl.read_csv(path, separator="\t", infer_schema_length=None)
            genes |= set(df[gene_col])
        elif path.suffix == ".csv":
            df = pl.read_csv(path, infer_schema_length=None)
            genes |= set(df[gene_col])
        else:
            with open(path, "r") as f:
                genes |= set(f.read().splitlines())
    _ = get_therapy_df(
        genes,
        cache=Path(cache),
        not_in_db_file=not_in_db,
        filter_confident=filter_confident,
        lookup=True,
    )


@click.command()
@click.option(
    "-i",
    "--cnv-input",
    required=True,
    help="""CNV file input, in the format <caller>:<file_path>.
    Can be specified multiple times to merge CNV calls from multiple
    callers
    """,
    default=None,
    multiple=True,
)
@click.option(
    "-c",
    "--concordant-out",
    required=True,
    help="Output file name for concordant calls",
    default=None,
)
@click.option(
    "-d",
    "--discordant-out",
    required=True,
    help="Output file name for discordant calls",
    default=None,
)
def merge_cnv_calls(
    cnv_input: list[str], concordant_out: str, discordant_out: str
) -> None:
    import polars as pl

    from sarek_plus.report.format import (
        AscatResult,
        CNVKitResult,
        CopyNumberResult,
        FreecResult,
    )

    caller_mapping: dict[str, Callable[[], CopyNumberResult]] = {
        "ascat": AscatResult,
        "cnvkit": CNVKitResult,
        "freec": FreecResult,
    }
    dfs: list[pl.DataFrame] = []
    for spec in cnv_input:
        caller, file = spec.split(":")
        if caller not in caller_mapping:
            raise ValueError(f"caller {caller} not recognized")
        cls = caller_mapping[caller]
        df = cls().read_segmental(file)
        dfs.append(df)
    if len(dfs) == 1:
        dfs[0].write_csv(concordant_out)
        dfs[0].filter(False).write_csv(discordant_out)
    else:
        concordant, discordant = CopyNumberResult.combine(dfs, add_genes=True)
        concordant.write_csv(concordant_out)
        discordant.write_csv(discordant_out)


@click.command()
@click.option("-i", "--input", required=True, help="Input file")
@click.option("-o", "--output", required=True, help="Output")
@click.option("-c", "--caller", required=True, help="Tool producing CNV calls")
@click.option(
    "-t", "--tumor_sample", required=False, help="Name of tumor sample", default=""
)
def classify_cnv_format(caller: str, input: str, output: str, tumor_sample: str):
    import polars as pl

    def dup_or_del(x):
        return "DEL" if x < 2 else "DUP"

    if caller.lower() == "cnvkit":
        df = pl.read_csv(input, separator="\t", infer_schema_length=None).select(
            ["chromosome", "start", "end", "cn"]
        )
    elif caller.lower() == "dellycnv":
        from subprocess import run

        if not tumor_sample:
            raise ValueError("Tumor sample name must be provided")
        qstring = f"bcftools query -f '%CHROM\t%POS0\t%END\t[%CN]\n' -s {tumor_sample} {input} > tmp.bed"
        run(
            qstring,
            shell=True,
        )
        df = pl.read_csv(
            "tmp.bed", separator="\t", new_columns=["chromosome", "start", "end", "cn"]
        )
        run("rm tmp.bed", shell=True)
    else:
        raise ValueError("Unknown caller provided")

    df.filter(pl.col("cn") != 2).with_columns(
        type=pl.col("cn").map_elements(dup_or_del, return_dtype=pl.String)
    ).drop("cn").write_csv(output, separator="\t", include_header=False)
