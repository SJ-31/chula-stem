#!/usr/bin/env python3
import io
import re
from functools import reduce
from subprocess import CompletedProcess, run

import click
import polars as pl
import polars.selectors as cs
from chula_stem.utils import empty_string2null


def in_vcf_header(vcf: str, string: str) -> bool:
    return run(f"bcftools head {vcf} | grep '{string}'", shell=True).returncode == 0


def _format_vep_vcf(
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
    df = _format_vep_vcf(
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


def _extract_format(vcf: str, samples: list[str] | None) -> pl.DataFrame:
    all_samples = (
        run(f"bcftools query -l '{vcf}'", shell=True, capture_output=True)
        .stdout.decode()
        .splitlines()
    )
    samples = list(set(samples) & set(all_samples)) if samples else all_samples
    if not samples:
        raise ValueError(
            f"None of the specified samples are in the VCF. Samples present: {all_samples}"
        )
    pos_df: pl.DataFrame = pl.read_csv(
        run(
            f"bcftools query -f '%CHROM\t%POS\t%ALT\t%REF' '{vcf}'",
            shell=True,
            capture_output=True,
        ).stdout,
        has_header=False,
        new_columns=["CHROM", "POS", "ALT", "REF"],
        separator="\t",
    )
    dfs = []
    for sample in samples:
        tmp = []
        for fline in (
            run(
                f"bcftools query -f '%FORMAT' '{vcf}' -s {sample}",
                shell=True,
                capture_output=True,
            )
            .stdout.decode()
            .splitlines()
        ):
            fields, vals = fline.split("\t")
            tmp.append({f: v for f, v in zip(fields.split(":"), vals.split(":"))})
        val_df = pl.from_dicts(tmp, strict=False, infer_schema_length=None)
        val_df = pl.concat([pos_df, val_df], how="horizontal").with_columns(
            pl.lit(sample).alias("sample")
        )
        dfs.append(val_df)
    return pl.concat(dfs, how="diagonal_relaxed")


# "^##FORMAT.*ID=\(.*?\),.*Description="\(.*\)">"


CALLER_FORMAT_FIELDS: dict[str, list] = {
    # Last updated [2026-11-10 Tue]
    "deepvariant": [
        "GT",  # Genotype
        "GQ",  # Conditional genotype quality
        "DP",  # Read depth
        "MIN_DP",  # Minimum DP observed within the GVCF block.
        "AD",  # Read depth for each allele
        "VAF",  # Variant allele fractions.
        "PL",  # Phred-scaled genotype likelihoods rounded to the closest integer
        "PS",  # Phase set
        "MF",  # Methylation fraction for each of the reference and alternate allele
        "MD",  # Methylation depth for each of the reference and alternate allele
        "MT",  # Methylation type: 0/0=Unmethylated, 0/1=Heterozygous, 1/1=Methylated
        "MI",  # Allele-specific methylation score: p-value for Wilcoxon Rank-Sum test based on the observed difference in methylation between haplotypes.
        "MED_DP",  # Median DP observed within the GVCF block rounded to the nearest integer.
    ],
    "manta_somatic": [
        "PR",  # Spanning paired-read support for the ref and alt alleles in the order listed
        "SR",  # Split reads for the ref and alt alleles in the order listed, for reads where P(allele|read)>0.999
    ],
    "haplotypecaller": [
        "AD",  # Allelic depths for the ref and alt alleles in the order listed
        "DP",  # Approximate read depth (reads with MQ=255 or with bad mates are filtered)
        "GQ",  # Genotype Quality
        "GT",  # Genotype
        "MIN_DP",  # Minimum DP observed within the GVCF block
        "PGT",  # Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles
        "PID",  # Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group
        "PL",  # Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification
        "PS",  # Phasing set (typically the position of the first variant in the set)
        "RGQ",  # Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)
        "SB",  # Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.
    ],
    "manta_diploid": [
        "GT",  # Genotype
        "FT",  # Sample filter, 'PASS' indicates that all filters have passed for this sample
        "GQ",  # Genotype Quality
        "PL",  # Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification
        "PR",  # Spanning paired-read support for the ref and alt alleles in the order listed
        "SR",  # Split reads for the ref and alt alleles in the order listed, for reads where P(allele|read)>0.999
    ],
    "mutect2": [
        "AD",  # Allelic depths for the ref and alt alleles in the order listed
        "AF",  # Allele fractions of alternate alleles in the tumor
        "DP",  # Approximate read depth (reads with MQ=255 or with bad mates are filtered)
        "F1R2",  # Count of reads in F1R2 pair orientation supporting each allele
        "F2R1",  # Count of reads in F2R1 pair orientation supporting each allele
        "FAD",  # Count of fragments supporting each allele.
        "GQ",  # Genotype Quality
        "GT",  # Genotype
        "PGT",  # Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles
        "PID",  # Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group
        "PL",  # Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification
        "PS",  # Phasing set (typically the position of the first variant in the set)
        "SB",  # Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.
    ],
    "freebayes": [
        "GT",  # Genotype
        "GQ",  # Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype
        "GL",  # Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy
        "DP",  # Read Depth
        "AD",  # Number of observation for each allele
        "RO",  # Reference allele observation count
        "QR",  # Sum of quality of the reference observations
        "AO",  # Alternate allele observation count
        "QA",  # Sum of quality of the alternate observations
        "MIN_DP",  # Minimum depth in gVCF output block.
    ],
    "strelka_snv": [  # Same for indels
        "DP",  # Read depth for tier1 (used+filtered)
        "FDP",  # Number of basecalls filtered from original read depth for tier1
        "SDP",  # Number of reads with deletions spanning this site at tier1
        "SUBDP",  # Number of reads below tier1 mapping quality threshold aligned across this site
        "AU",  # Number of 'A' alleles used in tiers 1,2
        "CU",  # Number of 'C' alleles used in tiers 1,2
        "GU",  # Number of 'G' alleles used in tiers 1,2
        "TU",  # Number of 'T' alleles used in tiers 1,2
    ],
    "muse": [
        "GT",  # Genotype
        "DP",  # Read depth at this position in the sample
        "AD",  # Depth of reads supporting alleles 0/1/2/3...
        "BQ",  # Average base quality for reads supporting alleles
        "SS",  # Variant status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post-transcriptional modification,5=unknown
    ],
    "deepsomatic": [
        "GT",  # Genotype
        "GQ",  # Conditional genotype quality
        "DP",  # Read depth
        "MIN_DP",  # Minimum DP observed within the GVCF block.
        "AD",  # Read depth for each allele
        "VAF",  # Variant allele fractions.
        "PL",  # Phred-scaled genotype likelihoods rounded to the closest integer
        "PS",  # Phase set
        "MF",  # Methylation fraction for each of the reference and alternate allele
        "MD",  # Methylation depth for each of the reference and alternate allele
        "MT",  # Methylation type: 0/0=Unmethylated, 0/1=Heterozygous, 1/1=Methylated
        "MI",  # Allele-specific methylation score: p-value for Wilcoxon Rank-Sum test based on the observed difference in methylation between haplotypes.
        "MED_DP",  # Median DP observed within the GVCF block rounded to the nearest integer.
        "NDP",  # Number of reads in the normal sample.
        "NAD",  # Read depth in the normal sample for alleles reported in the tumor sample
        "NAF",  # VAF of ALT alleles in the normal sample.
    ],
    "tiddit": [
        "CN",  # Copy number genotype for imprecise events
        "COV",  # Coverage (at A,B, and between)
        "DR",  # Number of paired-ends that supporting the reference allele (breakpoint A, and B)
        "DV",  # Number of paired-ends that support the event
        "GT",  # Genotype
        "LQ",  # Fraction of low quality reads
        "RR",  # Number of reads supporting the reference allele (breakpoint A, and B)
        "RV",  # Number of split reads that support the event
    ],
}


def merge_extracted(
    dfs: list[pl.DataFrame], callers: list[str] | None = None
) -> pl.DataFrame:
    callers = callers or [""] * len(dfs)
    tmp = []
    id_cols = ["CHROM", "POS", "ALT", "REF", "sample"]
    cols_merge = {("GQ", "mean"): set(), ("DP", "mean"): set(), ("GT", "unique"): set()}
    for df, caller in zip(dfs, callers):
        caller = f"{caller}_" if caller else ""
        df: pl.DataFrame = df.drop("caller") if "caller" in df.columns else df
        for (col, action), col_set in cols_merge.items():
            if col in df.columns:
                col_set.add(f"{caller}{col}")
                if action == "mean":
                    df = df.cast({col: pl.Float64})
        tmp.append(df.rename(lambda c: c if c in id_cols else f"{caller}{c}"))
    merged: pl.DataFrame = reduce(
        lambda x, y: x.join(y, how="full", on=id_cols)
        .with_columns(pl.coalesce(col, f"{col}_right") for col in id_cols)
        .drop(cs.ends_with("_right")),
        tmp,
    )
    for (col, action), col_set in cols_merge.items():
        if not col_set:
            continue
        name = f"{col}_{action}"
        if action == "unique":
            col_add = pl.concat_list(col_set).list.unique().list.join(",")
        elif action == "mean":
            col_add = merged.select(col_set).mean_horizontal()
        merged = merged.with_columns(col_add.alias(name))
    return merged


@click.command()
@click.option(
    "-i",
    "--input",
    required=True,
    help="""Input VCF file, or a string <CALLER>:<VCF_FILENAME>.
    Can be specified multiple times for files from separate callers.
    The format fields will be merged appropriately
    """,
    multiple=True,
    default=None,
)
@click.option(
    "-s",
    "--sample",
    required=False,
    help="Name of sample to extract from file. Defaults to all samples if not provided",
    default=None,
    multiple=True,
)
@click.option(
    "-r",
    "--samples_rename",
    required=False,
    help="Samples to rename, written as OLD_NAME=NEWNAME pairs",
    default=None,
    multiple=True,
)
@click.option("-o", "--output", required=True, help="Output TSV file", default=None)
@click.option("-c", "--caller", required=False, help="Name of caller", default=None)
def extract_vcf_format(
    input: tuple[str],
    output: str,
    caller: str,
    sample,
    samples_rename: tuple[str] | None = None,
) -> None:
    to_rename: dict = (
        {o[0]: o[1] for o in (s.split("=") for s in samples_rename)}
        if samples_rename
        else {}
    )
    if len(input) == 1:
        df: pl.DataFrame = _extract_format(input[0], samples=sample).with_columns(
            pl.col("sample").replace(to_rename)
        )
        if caller:
            df = df.with_columns(pl.lit(caller).alias("caller"))
    else:
        tmp, callers = [], []
        for cur in input:
            cur_caller = None
            if ":" in cur:
                cur_caller, cur = cur.split(":")
            extracted = _extract_format(cur, samples=sample).with_columns(
                pl.col("sample").replace(to_rename)
            )
            if cur_caller:
                extracted = extracted.with_columns(pl.lit(cur_caller).alias("caller"))
                callers.append(cur_caller)
            tmp.append(extracted)
        df = merge_extracted(tmp, callers or None)
    df.write_csv(output, separator="\t")
