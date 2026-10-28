#!/usr/bin/env python3
import io
import re
from subprocess import CompletedProcess, run

import click
import polars as pl
import polars.selectors as cs
from chula_stem.utils import empty_string2null


def in_vcf_header(vcf: str, string: str) -> bool:
    return run(f"bcftools head {vcf} | grep '{string}'", shell=True).returncode == 0


def format_vep_vcf_internal(
    input: str,
    tumor_sample: str,
    normal_sample: str = "",
    vep_info_field: str = "ANN",
    tool_source_tag: str | None = "TOOL_SOURCE",
    vaf_tag: str = "AF",
    ad_tag: str = "AD",
    variant_class: str = "small",
) -> pl.DataFrame:
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
    if tool_source_tag is None:
        vcf_cols.remove(rf"%INFO/{tool_source_tag}")
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

    if tool_source_tag is not None:
        first_cols = ["Loc", "Ref", "Alt", tool_source_tag, "FILTER"]
    else:
        first_cols = ["Loc", "Ref", "Alt", "FILTER"]
    wanted_cols = first_cols + vaf_ad_cols + vep_columns
    if variant_class == "sv":
        wanted_cols.append("SVTYPE")
    df = df.select(wanted_cols)
    return df


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
def main(
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
    df = format_vep_vcf_internal(
        input=input,
        tumor_sample=tumor_sample,
        normal_sample=normal_sample,
        vep_info_field=vep_info_field,
        tool_source_tag=tool_source_tag,
        vaf_tag=vaf_tag,
        ad_tag=ad_tag,
        variant_class=variant_class,
    )
    df.write_csv(output, separator="\t", null_value="NA")
