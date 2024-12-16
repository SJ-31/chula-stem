#!/usr/bin/env ipython
# import chula_stem as
import pytest
from chula_stem.utils import _format_vep_vcf

vepdir = "/home/shannc/Bio_SDD/chula-stem/tests/vep"


@pytest.mark.skip(reason="Done")
def test_format_vep():
    sv = "/home/shannc/Bio_SDD/chula-stem/tests/vep/7-patient_10-VEP_SV.vcf.gz"
    small = (
        "/home/shannc/Bio_SDD/chula-stem/tests/vep/7-patient_10-VEP_small.vcf.gz"
    )
    _format_vep_vcf(
        sv,
        f"{vepdir}/sv.tsv",
        tumor_sample="10_cancer",
        normal_sample="10_blood",
        tool_source_tag="SOURCE",
    )
    _format_vep_vcf(
        small,
        f"{vepdir}/small.tsv",
        tumor_sample="10_cancer",
        normal_sample="10_blood",
        tool_source_tag="SOURCE",
    )
    _format_vep_vcf(
        small,
        f"{vepdir}/small_no.tsv",
        tumor_sample="10_cancer",
        tool_source_tag="SOURCE",
    )


##
import polars as pl

msidir = "/home/shannc/Bio_SDD/chula-stem/tests/msisensor"
unstable_path = f"{msidir}/4-null-Msisensor_unstable.tsv"
overlapping_path = f"{msidir}/unique.bed"

unstable: pl.DataFrame = (
    pl.read_csv(unstable_path, separator="\t", infer_schema_length=None)
    .with_columns(
        end=pl.col("location")
        + (pl.col("repeat_unit_bases").str.len_chars() * pl.col("repeat_times"))
    )
    .rename({"location": "start"})
)
bed_cols: list = ["chromosome", "start", "end", "gene_name"]
overlapping: pl.DataFrame = pl.read_csv(
    overlapping_path,
    separator="\t",
    infer_schema_length=None,
    new_columns=bed_cols,
)
unstable.join_where(
    overlapping,
    pl.col("start") >= pl.col("start_right"),
    pl.col("end") <= pl.col("end_right"),
)
