#!/usr/bin/env ipython
# import chula_stem as

import pytest
from chula_stem.utils import _format_vep_vcf
from rpy2.robjects.packages import STAP, importr

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


import pandas as pd
import polars as pl
from rpy2.robjects import pandas2ri


def pd2r(df: pl.DataFrame | pd.DataFrame):
    if isinstance(df, pl.DataFrame):
        df = df.to_pandas()
    with (ro.default_converter + pandas2ri.converter).context():
        return ro.conversion.get_conversion().py2rpy(df)


## want to compare your called mutations against the previous ones for PDAC
import rpy2.robjects as ro
from chula_stem.utils import r2pd
from rpy2.robjects.packages import importr

classify_cnv = (
    "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/8-null-ClassifyCNV.tsv"
)
ref = "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/aggregated_cnv.tsv"

results = (
    pl.read_csv(
        classify_cnv, separator="\t", null_values="NA", infer_schema_length=None
    )
    .filter(~pl.any_horizontal(pl.col(["Start", "End"]).str.contains("[a-zA-Z]")))
    .cast({"Start": pl.Int64, "End": pl.Int64})
)
ag = pl.read_csv(ref, separator="\t", null_values="NA", infer_schema_length=None)
dplyr = importr("dplyr")

data_r = pd2r(results)
ref_r = pd2r(ag)
# joined = dplyr.inner_join(data_r, ref_r, by=dplyr.join_by("chr", overlaps()))
# inner_join(cnv, ref, by = join_by(x$Chromosome == y$chr, overlaps(x$Start, x$End, y$start, y$stop)))

##
# TODO: see if you can use this
src = """
library(tidyverse)

overlapping_join <- function(
    x, y, x_on, y_on, x_start = "start", x_end = "end", y_start = "start",
    y_end = "end") {
  x_rename <- c(tempX = x_on, x_start = x_start, x_end = x_end)
  y_rename <- c(tempY = y_on, y_start = y_start, y_end = y_end)
  rename_back <- c(names(x_rename), names(y_rename)) |>
    `names<-`(c(x_rename, y_rename))
  x <- rename(x, all_of(x_rename))
  y <- rename(y, all_of(y_rename))
  inner_join(x, y, by = join_by(
    x$tempX == y$tempY,
    overlaps(x_start, x_end, y_start, y_end)
  )) |> rename(any_of(rename_back))
}

"""
