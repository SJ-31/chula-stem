#!/usr/bin/env ipython

import ast
import itertools
import re
from functools import reduce
from pathlib import Path
from subprocess import CompletedProcess, run
from tempfile import TemporaryFile
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


@click.command()
@click.option("-n", "--normal", required=True)
@click.option("-t", "--tumor", required=True)
@click.option("-s", "--snps", help="Path to snps file", required=True)
@click.option(
    "-r", "--chr_files", help="Path to directory with chromosomes", required=True
)
@click.option("-l", "--chr_len_file", required=True)
@click.option("-o", "--output", required=True)
@click.option("-b", "--breakpoint_threshold", default=1.2)
@click.option("-m", "--minimal_coverage_per_position", default=0)
@click.option("-p", "--ploidy", default="2,3")
@click.option(
    "-r",
    "--read_type",
    default="paired_end",
    type=click.Choice(["paired_end", "mate_pair", "single"]),
)
@click.option("-w", "--window", default=0)
@click.option("-d", "--outputdir", default=".")
@click.option("-c", "--contamination", default=0)
@click.option("-e", "--sex", default="XY")
@click.option("-y", "--breakpoint_type", default=2)
@click.option("-x", "--max_threads", default=1)
@click.option("--noisy_data", default=True, type=bool, is_flag=True)
@click.option(
    "-f",
    "--file_format",
    default="BAM",
    help="input format of tumor and normal files",
)
@click.option("--print_na", default=False, type=bool, is_flag=True)
@click.option("-i", "--intervals", help="Exome target intervals")
def make_freec_config(
    tumor: str,
    normal: str,
    chr_len_file: str,
    chr_files: str,
    max_threads: int,
    minimal_coverage_per_position,
    window,
    file_format,
    ploidy,
    outputdir,
    contamination,
    breakpoint_threshold,
    breakpoint_type,
    sex,
    snps,
    noisy_data,
    print_na,
    read_type,
    intervals,
    output,
):
    import toml

    template: dict = {
        "general": {
            "chrLenFile": chr_len_file,
            "window": window,
            "ploidy": ploidy,
            "outputDir": outputdir,
            "sex": sex,
            "breakPointType": breakpoint_type,
            "contamination": contamination,
            "chrFiles": chr_files,
            "maxThreads": max_threads,
            "breakPointThreshold": breakpoint_threshold,
        },
        "sample": {"mateFile": normal, "inputFormat": file_format},
        "control": {"mateFile": tumor, "inputFormat": file_format},
        "BAF": {
            "SNPfile": snps,
            "minimalCoveragePerPosition": minimal_coverage_per_position,
        },
    }
    if read_type == "paired_end":
        template["sample"]["mateOrientation"] = "FR"
        template["control"]["mateOrientation"] = "FR"
    elif read_type == "mate_pair":
        template["sample"]["mateOrientation"] = "RF"
        template["control"]["mateOrientation"] = "RF"
    elif read_type == "single":
        template["sample"]["mateOrientation"] = 0
        template["control"]["mateOrientation"] = 0
    if noisy_data:
        template["general"]["noisyData"] = "TRUE"
    if not print_na:
        template["general"]["printNA"] = "false"
    if intervals:
        template["target"] = {"captureRegions": intervals}
        template["general"]["window"] = 0
    with TemporaryFile("r+") as w:
        toml.dump(template, w)
        w.seek(0)
        lines = w.readlines()
    with open(output, "w") as w:
        modified: list[str] = [l.replace('"', "").replace("'", "") for l in lines]
        w.write("".join(modified))


def in_vcf_header(vcf: str, string: str) -> bool:
    from subprocess import run

    return run(f"bcftools head {vcf} | grep '{string}'", shell=True).returncode == 0


@click.command()
@click.option("-i", "--input", required=True, help="Input vcf file")
@click.option("-o", "--output", required=True, help="Output tsv")
@click.option(
    "-t",
    "--tool_source_tag",
    required=False,
    help="Info tag designating caller of origin",
)
@click.option(
    "-n",
    "--normal_sample",
    required=False,
    default="",
    help="Name of the normal sample in the vcf file",
)
@click.option(
    "-r",
    "--tumor_sample",
    required=True,
    help="Name of the tumor sample in the vcf file",
)
@click.option(
    "-c",
    "--variant_class",
    required=False,
    help="Type of variants in file",
    default="small",
    type=click.Choice(["small", "sv"]),
)
@click.option(
    "-v",
    "--vep_info_field",
    required=True,
    help="Info field in vcf containing VEP annotations",
)
@click.option("--vaf_tag", required=False, help="Name of VAF tag", default="AF")
@click.option("--ad_tag", required=False, help="Name of AD tag", default="AD")
def format_vep_vcf(
    input: str,
    output: str,
    tumor_sample: str,
    normal_sample: str = "",
    vep_info_field: str = "ANN",
    tool_source_tag: str = "TOOL_SOURCE",
    vaf_tag: str = "AF",
    ad_tag: str = "AD",
    variant_class: str = "small",
):
    import io
    import re

    import polars as pl
    import polars.selectors as cs

    proc: CompletedProcess = run(
        f"bcftools head {input} | grep 'INFO=<ID={vep_info_field}'",
        shell=True,
        capture_output=True,
        check=True,
    )
    vep_columns: list = re.findall('.*Format: (.*)">', proc.stdout.decode())[0].split(
        "|"
    )
    vaf: str = "[%AF\t]" if in_vcf_header(input, f"##FORMAT=<ID={vaf_tag}") else ""
    ad: str = "[%AD\t]" if in_vcf_header(input, f"##FORMAT=<ID={ad_tag}") else ""
    vcf_cols: list = [
        "%CHROM:%POS",
        "%REF",
        "%ALT",
        vaf,
        ad,
        rf"%INFO/{tool_source_tag}",
        "%FILTER",
        rf"%INFO/{vep_info_field}",
    ]
    vaf_ad_cols: list = []
    if vaf:
        vaf_ad_cols.append("VAF")
    if vaf and normal_sample:
        vaf_ad_cols.append("VAF_normal")
    if ad:
        vaf_ad_cols.append("Alt_depth")
    if ad and normal_sample:
        vaf_ad_cols.append("Alt_depth_normal")

    sample_flag: str = (
        f"{tumor_sample},{normal_sample}" if normal_sample else tumor_sample
    )
    columns = ["Loc", "Ref", "Alt"] + vaf_ad_cols + [tool_source_tag, "FILTER", "ANN"]
    if variant_class == "sv":
        columns.insert(3, "SVTYPE")
        vcf_cols.insert(3, "%INFO/SVTYPE")
    qstring = "\t".join(list(filter(lambda x: x, vcf_cols))).replace("]\t", "]")
    runstr = rf"bcftools query -f '{qstring}' -s '{sample_flag}' {input}"
    proc2: CompletedProcess = run(runstr, shell=True, capture_output=True, check=True)
    df = (
        pl.read_csv(
            io.StringIO(proc2.stdout.decode()),
            separator="\t",
            new_columns=columns,
            null_values=".",
            infer_schema_length=None,
        )
        .with_columns(pl.col("ANN").str.split(","))
        .explode("ANN")
        .with_columns(pl.col("ANN").str.split("|").list.to_struct(fields=vep_columns))
        .unnest("ANN")
        .pipe(empty_string2null)
    )
    if not vaf:
        df = df.with_columns(VAF=pl.lit(None))
    if not ad:
        df = df.with_columns(Alt_depth=pl.lit(None))
    else:
        alt_cols = list(filter(lambda x: "Alt_depth" in x, vaf_ad_cols))
        ad_split = [
            pl.col(x).str.split(",").list.get(1, null_on_oob=True) for x in alt_cols
        ]
        df = df.with_columns(ad_split).cast({a: pl.Int32 for a in alt_cols})
        if pl.String in df.select(cs.starts_with("VAF")).dtypes:
            replace_dots = cs.starts_with("VAF").str.replace(",.", "", literal=True)
            df = df.with_columns(replace_dots)

    wanted_cols = (
        ["Loc", "Ref", "Alt", tool_source_tag, "FILTER"] + vaf_ad_cols + vep_columns
    )
    if variant_class == "sv":
        wanted_cols.append("SVTYPE")
    df = df.select(wanted_cols)
    df.write_csv(output, separator="\t", null_value="NA")


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


def read_existing[T](filename: Path, expr: Callable[[Path], T], read_fn) -> T:
    if filename.exists():
        return read_fn(filename)
    else:
        return expr(filename)
