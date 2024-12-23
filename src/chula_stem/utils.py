#!/usr/bin/env ipython


from subprocess import CompletedProcess, run
from tempfile import TemporaryFile

import click
import pandas as pd
import polars as pl


# rx method is equivalent to [], rx2 is [[]]
def empty_string2null(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(
        pl.when(pl.col(pl.String).str.len_chars() == 0)
        .then(None)
        .otherwise(pl.col(pl.String))
        .name.keep()
    )


def r2pd(robject) -> pd.DataFrame:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    with (ro.default_converter + pandas2ri.converter).context():
        obj = ro.conversion.get_conversion().rpy2py(robject)
        if isinstance(obj, tuple) and isinstance(obj[0], pd.DataFrame):
            return obj[0]
        return obj


@click.command()
@click.option("-i", "--input", required=True, help="Input file")
@click.option("-o", "--output", required=True, help="Output")
@click.option("-c", "--caller", required=True, help="Tool producing CNV calls")
def classify_cnv_format(caller: str, input: str, output: str):
    _classify_cnv(caller, input, output)


def read_facets_rds(rds_path: str) -> pl.DataFrame:
    base = importr("base")
    df: pd.DataFrame = r2pd(base.readRDS(rds_path).rx2("segs"))
    try:
        df = df.astype({"chrom": "int32"})
    except:
        pass
    return pl.from_pandas(df, schema_overrides={"chrom": pl.String})


def _classify_cnv(caller: str, input: str, output: str, tumor_sample: str = ""):
    import polars as pl

    def dup_or_del(x):
        return "DEL" if x < 2 else "DUP"

    if caller.lower() == "cnvkit":
        df = pl.read_csv(input, separator="\t").select(
            ["chromosome", "start", "end", "cn"]
        )
    elif caller.lower() == "facets":
        base = importr("base")
        df: pl.DataFrame = (
            pl.from_pandas(r2pd(base.readRDS(input).rx2("segs")))
            .select(["chrom", "start", "end", "tcn.em"])
            .rename({"tcn.em": "cn"})
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
            "tmp.bed",
            separator="\t",
            new_columns=["chromosome", "start", "end", "cn"],
            dtypes=[pl.String, pl.Int64, pl.Int64, pl.Int64],
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
    "-v",
    "--vep_info_field",
    required=True,
    help="Info field in vcf containing VEP annotations",
)
def format_vep_vcf(
    input, output, tumor_sample, normal_sample, vep_info_field, tool_source_tag
):
    _format_vep_vcf(
        input,
        output,
        tumor_sample,
        normal_sample,
        vep_info_field,
        tool_source_tag,
    )


def in_vcf_header(vcf: str, string: str) -> bool:
    from subprocess import run

    return (
        run(f"bcftools head {vcf} | grep '{string}'", shell=True).returncode == 0
    )


def _format_vep_vcf(
    input: str,
    output: str,
    tumor_sample: str,
    normal_sample: str = "",
    vep_info_field: str = "ANN",
    tool_source_tag: str = "TOOL_SOURCE",
    vaf_tag: str = "AF",
    ad_tag: str = "AD",
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
    vep_columns: list = re.findall('.*Format: (.*)">', proc.stdout.decode())[
        0
    ].split("|")
    vaf: str = (
        "[%AF\t]" if in_vcf_header(input, f"##FORMAT=<ID={vaf_tag}") else ""
    )
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
    qstring = "\t".join(list(filter(lambda x: x, vcf_cols))).replace("]\t", "]")
    runstr = rf"bcftools query -f '{qstring}' -s '{sample_flag}' {input}"
    proc2: CompletedProcess = run(
        runstr, shell=True, capture_output=True, check=True
    )
    columns = (
        ["Loc", "Ref", "Alt"] + vaf_ad_cols + [tool_source_tag, "FILTER", "ANN"]
    )
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
        .with_columns(
            pl.col("ANN").str.split("|").list.to_struct(fields=vep_columns)
        )
        .unnest("ANN")
        .pipe(empty_string2null)
    )
    if not vaf:
        df = df.with_columns(VAF=pl.lit(None))
    if not ad:
        df = df.with_columns(Alt_depth=pl.lit(None))
    else:
        replace_dots = cs.starts_with("VAF").str.replace(",.", "", literal=True)
        alt_cols = list(filter(lambda x: "Alt_depth" in x, vaf_ad_cols))
        ad_split = [
            pl.col(x).str.split(",").list.get(1, null_on_oob=True)
            for x in alt_cols
        ]
        df = (
            df.with_columns(ad_split)
            .cast({a: pl.Int32 for a in alt_cols})
            .with_columns(replace_dots)
        )
    df = df.select(
        ["Loc", "Ref", "Alt", tool_source_tag, "FILTER"]
        + vaf_ad_cols
        + vep_columns
    )
    df.write_csv(output, separator="\t", null_value="NA")
