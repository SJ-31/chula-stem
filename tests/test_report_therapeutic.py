#!/usr/bin/env ipython

from collections import Counter
from os import replace

import polars as pl
import polars.selectors as cs
from chula_stem.report.format import dbvar_link, get_clingen_link
from chula_stem.utils import read_facets_rds
from chula_stem.report import ReportElement, ResultsReport
from chula_stem.report.spec import URL, Rename, Widths
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm, inch
from reportlab.pdfgen.canvas import Canvas

## Working on formatting the therapeutics table
table_spec = {"header_pos": (A4[0] - 5 * inch, A4[1] - inch)}
numeric_style: ParagraphStyle = ParagraphStyle("nums", fontSize=11)
text_style: ParagraphStyle = ParagraphStyle("data", fontSize=10)

# TODO: add a style for underlining links and highlighting them
cell_pstyles: dict = {
    1: numeric_style,
    2: numeric_style,
    None: text_style,
}
PUBCHEM_URL: str = "https://pubchem.ncbi.nlm.nih.gov/compound"

cell_styles = style_cells((0, 1), background=colors.lightcyan, valign="TOP")
header_styles = style_cells(
    (0, 0),
    8,
    1,
    textcolor=colors.red,
    underline=(3, colors.black),
    background=colors.lightgrey,
)

ttable_styles: dict = {
    "cell_pstyles": cell_pstyles,
    "header_pstyles": ParagraphStyle("cols", fontSize=9),
    "cell_styles": cell_styles,
    "header_styles": header_styles,
    "col_widths": list(Widths.therapy.values()),
}


TUMOR_KEYWORDS = ["cancer", "leukemia", "carcinoma", "lymphoma"]

# Will have
# May have to have tables in tables to show the therapy info properly

# For small and structural variants
# TODO: will also need to do this for mutational signatures.
# CNV data is a bonus, but probably not possible
#

#
# others = df.filter(pl.col("db") != "pandrugs2")
# * Treatment
# # pubchemId
# * Evidence type
#   * Can be mutational signature or variant
# * Relevant evidence in sample (gene names or mutational signature ids)
# * Use in other cancers
# *  Study
# # Db link
#

## Repetitive elements
from chula_stem.utils import add_loc
from chula_stem.report.spec import Rename

file = "/home/shannc/Bio_SDD/chula-stem/tests/msisensor/4-null-CR.tsv"
df = pl.read_csv(file, separator="\t", null_values="NA", infer_schema_length=None)
df = (
    add_loc(df, start_col="Start", end_col="End")
    .with_columns(
        pl.col("ClinGen_report").map_elements(get_clingen_link, return_dtype=pl.String)
    )
    .with_columns(
        pl.struct(s="source", acc="accession")
        .map_elements(
            lambda x: (
                dbvar_link(x["acc"]) if x["acc"] != "NA" and x["acc"] else x["s"]
            ),
            return_dtype=pl.String,
        )
        .alias("source")
    )
).rename(Rename.repeat)


small = pl.read_parquet(
        "/home/shannc/Bio_SDD/chula-stem/tests/vep_format_sv/small_all.parquet"
)
rel = pl.read_parquet(
    "/home/shannc/Bio_SDD/chula-stem/tests/vep_format_sv/small_relevant.parquet"
)
sv = pl.read_parquet(
    "/home/shannc/Bio_SDD/chula-stem/tests/vep_format_sv/sv_all.parquet"
)
