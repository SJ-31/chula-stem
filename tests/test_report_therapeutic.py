#!/usr/bin/env ipython

from collections import Counter
from os import replace

import polars as pl
import polars.selectors as cs
from chula_stem.utils import read_facets_rds
from chula_stem.report import ReportElement, ResultsReport, style_cells, add_link
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
all = "/home/shannc/Bio_SDD/chula-stem/tests/report_tmp/small_all.parquet"

# For small and structural variants
# TODO: will also need to do this for mutational signatures.
# CNV data is a bonus, but probably not possible
#


db_link_fn = lambda x, y: add_link(x, y, underline="yes", underlinecolor="#5e81ac")
df = (
    (
        (
            pl.read_parquet(all)
            .select(pl.col(["Gene", "disease", "source", "therapies", "db", "db_link"]))
            .explode("disease")
        )
        .with_columns(
            pl.col("disease").str.to_lowercase().str.replace_all("_", " "),
            pl.struct(["db", "db_link"])
            .map_elements(
                lambda x: db_link_fn(x["db"], x["db_link"]),
                return_dtype=pl.String,
            )
            .alias("db_link"),
        )
        .filter(  # Make sure that the therapies reported are relevant to cancer
            (pl.col("db") == "pandrugs2")
            | pl.col("disease").str.contains_any(TUMOR_KEYWORDS)
        )
        .with_columns(
            pl.col("disease")
            .str.replace_many({"cancer": "", "clinical": ""})
            .str.replace("", "unspecified")
            .map_elements(str.capitalize, return_dtype=pl.String)
        )
    )
    .explode("therapies")
    .with_columns(
        pl.col("therapies")
        .str.split(":")
        .list.to_struct(fields=["therapies", "PubChemId"])
    )
    .unnest("therapies")
    .with_columns(
        pl.col("PubChemId").map_elements(
            lambda x: (
                add_link(x, f"{URL.pubchem}/{x}", underline="yes") if x != "NA" else x
            ),
            return_dtype=pl.String,
        ),
        pl.lit("Gene variation").alias("type"),
    )
    .filter(pl.col("therapies").str.contains("[A-Z-a-z]*"))
    .rename(Rename.therapy)
).select(list((Rename.therapy).values()))

ResultsReport.build_table(
    table_spec,
    ttable_styles,
    df,
    "Relevant therapies",
    "Relevant therapies",
    None,
    "/home/shannc/Bio_SDD/chula-stem/tests/therapy_table.pdf",
)

# "therapy" will be the unique column
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
