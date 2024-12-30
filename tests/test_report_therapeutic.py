#!/usr/bin/env ipython

from collections import Counter
import re
from os import replace

import polars as pl
import polars.selectors as cs
from chula_stem.report.format import (
    add_link,
    dbvar_link,
    get_clingen_link,
    sigprofiler_fmt,
)
from chula_stem.utils import read_facets_rds
from chula_stem.report import ReportElement, ResultsReport
from chula_stem.report.spec import URL, Rename, Widths, signature_style
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm, inch
from reportlab.pdfgen.canvas import Canvas

from resources.cosmic_signatures import COSMIC_URL

## Working on formatting the therapeutics table
table_spec = {"header_pos": (A4[0] - 5 * inch, A4[1] - inch)}
numeric_style: ParagraphStyle = ParagraphStyle("nums", fontSize=11)
text_style: ParagraphStyle = ParagraphStyle("data", fontSize=10)

# Will have
# May have to have tables in tables to show the therapy info properly

# For small and structural variants
# TODO: will also need to do this for mutational signatures.
# CNV data is a bonus, but probably not possible
#

#
# others = df.filter(pl.col("db") != "pandrugs2")
# - Treatment
# # pubchemId
# - Evidence type
#   * Can be mutational signature or variant
# - Relevant evidence in sample (gene names or mutational signature ids)
# - Use in other cancers
# -  Study
# # Db link
#
from pypdf import PdfWriter
from pypdf.annotations import FreeText
import pymupdf
from pymupdf import Document, Page


def add_page_numbers(page: int, total: int, writer: PdfWriter, offset: int = 1):
    # Rectangle coordinates are [x1, y1, x2, y2]
    # Where x1 and y1 are the coordinates of the bottom-left corner
    # And x2 and y2 are the coordinates of the top-right corner
    text = FreeText(
        text=f"Page {page + offset} of {total}",
        font="Courier",
        font_size="15pt",
        rect=(A4[0] - 3 * inch, A4[1] - inch, A4[0] - inch, A4[1]),
        bold=True,
        italic=True,
        background_color=None,
        border_color=None,
    )
    writer.add_annotation(page_number=page, annotation=text)


reportdir = "/home/shannc/Bio_SDD/chula-stem/tests/report_tmp"

tables = ["small", "sv", "cnv", "repeat"]
names = [
    "Small Variants",
    "Structural Variants",
    "Copy Number Variants",
    "Tandem Repeats",
]
# * Therapy formatter with signatures
from chula_stem.report.format import therapy_fmt, sigprofiler_fmt

data = pl.read_parquet(f"{reportdir}/therapy_data.parquet")
