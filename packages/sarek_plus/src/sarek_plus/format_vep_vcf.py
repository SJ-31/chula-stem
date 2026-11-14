#!/usr/bin/env python3
import io
import re
from functools import reduce
from io import StringIO
from pathlib import Path
from subprocess import CompletedProcess, run
from typing import Literal

try:
    from icecream import ic

    ic.configureOutput("dbg: ", includeContext=True)
except ImportError:  # Graceful fallback if IceCream isn't installed.
    ic = lambda *a: None if not a else (a[0] if len(a) == 1 else a)  # noqa
import click
import polars as pl
import polars.selectors as cs
from chula_stem.utils import empty_string2null


def in_vcf_header(vcf: str, string: str) -> bool:
    return run(f"bcftools head {vcf} | grep '{string}'", shell=True).returncode == 0


def _format_vep_vcf(
    input: str,
    vep_info_field: str = "ANN",
    fields: tuple[str] | None = None,
) -> pl.DataFrame:
    proc: CompletedProcess = run(
        f"bcftools head {input} | grep 'INFO=<ID={vep_info_field}'",
        shell=True,
        capture_output=True,
        check=True,
    )
    info_fields = [
        re.findall(".*ID=(.*?),", line)[0]
        for line in (
            run(
                f"bcftools head '{input}' | grep '^##INFO'",
                shell=True,
                capture_output=True,
                check=True,
            )
            .stdout.decode()
            .splitlines()
        )
    ]
    info_fields.remove(vep_info_field)
    fields_keep = list(set(info_fields) & set(fields)) if fields else []

    vep_columns: list = re.findall('.*Format: (.*)">', proc.stdout.decode())[0].split(
        "|"
    )
    vcf_cols: list = [
        "%CHROM",
        "%POS",
        "%REF",
        "%ALT",
        "%FILTER",
        rf"%INFO/{vep_info_field}",
    ] + [rf"%INFO/{field}" for field in fields_keep]
    columns = ["CHROM", "POS", "REF", "ALT"] + ["FILTER", vep_info_field] + fields_keep
    qstring = "\t".join(list(filter(lambda x: x, vcf_cols)))
    runstr = rf"bcftools query -f '{qstring}' {input}"
    proc2: CompletedProcess = run(runstr, shell=True, capture_output=True, check=True)
    df = (
        pl.read_csv(
            io.StringIO(proc2.stdout.decode()),
            separator="\t",
            new_columns=columns,
            null_values=".",
            infer_schema_length=None,
        )
        .with_columns(pl.col(vep_info_field).str.split(","))
        .explode(vep_info_field)
        .with_columns(
            pl.col(vep_info_field).str.split("|").list.to_struct(fields=vep_columns)
        )
        .unnest(vep_info_field)
        .pipe(empty_string2null)
    )
    first_cols = ["CHROM", "POS", "REF", "ALT"]
    wanted_cols = first_cols + vep_columns + ["FILTER"] + fields_keep
    df = df.select(wanted_cols)
    return df


@click.command()
@click.option("-i", "--input", required=True, help="Input vcf file")
@click.option("-o", "--output", required=True, help="Output tsv")
@click.option(
    "-f",
    "--fields",
    required=False,
    help="Additional INFO fields to retrieve",
    default=None,
    multiple=True,
)
@click.option(
    "-v",
    "--vep_info_field",
    required=True,
    help="Info field in vcf containing VEP annotations",
)
def main(
    input: str,
    output: str,
    vep_info_field: str = "ANN",
    fields: tuple[str] | None = None,
):
    df = _format_vep_vcf(
        input=input,
        vep_info_field=vep_info_field,
        fields=fields,
    )
    df.write_csv(output, separator="\t", null_value="NA")


def get_vcf_types(
    vcf: str, field: Literal["FORMAT", "INFO", "FILTER"] = "FORMAT"
) -> dict:
    pl_type_mapping = {
        "Integer": pl.Int64,
        "Number": pl.Float64,
        "Float": pl.Float64,
        "Flag": pl.Boolean,
        "Character": pl.String,
        "String": pl.String,
    }
    type_mapping = {
        k: pl_type_mapping[v] if n == "1" else pl.String
        for k, n, v in (
            re.findall("ID=(.*?),Number=(.*?),Type=(.*?),", line)[0]
            for line in run(f"bcftools head '{vcf}'", capture_output=True, shell=True)
            .stdout.decode()
            .splitlines()
            if line.startswith(f"##{field}")
        )
    }
    return type_mapping


def _extract_format(vcf: str, samples: list[str] | None) -> pl.DataFrame:
    assert Path(vcf).exists(), f"File `{vcf}` doesn't exist"
    all_samples = (
        run(f"bcftools query -l '{vcf}'", shell=True, capture_output=True)
        .stdout.decode()
        .splitlines()
    )
    samples = list(set(samples) & set(all_samples)) if samples else all_samples
    if not samples:
        raise ValueError(
            f"None of the specified samples are in the VCF `{vcf}`. Samples present: {all_samples}"
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
        null_values=["."],
    )
    format_types = get_vcf_types(vcf, "FORMAT")
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
        val_df = val_df.with_columns(cs.by_dtype(pl.String).replace(".", None)).cast(
            {k: v for k, v in format_types.items() if k in val_df.columns}
        )
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
    df.write_csv(output, separator="\t", null_value="NA")
