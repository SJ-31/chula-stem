from pathlib import Path
from typing import Callable

import chula_stem.report as rp
import polars as pl
from chula_stem.callset_qc import IMPACT_MAP
from chula_stem.report import VTABLE_COL_WIDTHS, Civic, PanDrugs2, TherapyDB

small_path = (
    "/home/shannc/Bio_SDD/chula-stem/tests/vep/7-patient_10-VEP_small_1.tsv"
)
sv_pat = "/home/shannc/Bio_SDD/chula-stem/tests/vep/7-patient_10-VEP_SV_1.tsv"

# Documents https://docs.reportlab.com/reportlab/userguide/ch5_platypus/#documents-and-templates
# Create a list of flowables and pass it to Doc.build()
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm, inch
from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus import (
    BaseDocTemplate,
    LongTable,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
)

STYLES = getSampleStyleSheet()
PWIDTH = A4[0]
PHEIGHT = A4[1]

# (0, 0) is the bottom left of the page by default

# Flowables
# - Flowable.drawOn(canvas, x, y)

# You can work in cm by converting doc.width (or height) with the `units` module like so:
# doc.width/units.cm

sample = "/home/shannc/Bio_SDD/chula-stem/tests/report_sample.tsv"
if not Path(sample).exists():
    R = rp.ResultsReport(
        "dummy",
        civic_cache="/home/shannc/Bio_SDD/chula-stem/tests/civic.json",
        pandrugs2_cache="/home/shannc/Bio_SDD/chula-stem/tests/pandrugs2.json",
        vep_small=small_path,
        vep_sv=sv_pat,
    )
    res = R.data["all"]["small"]
    res.write_csv(sample, separator="\t", null_value="NA")
else:
    res = pl.read_csv(sample, separator="\t", null_values="NA")


# Use a doc template to enable the table to be split across multiple pages
out = "/home/shannc/Bio_SDD/chula-stem/test.pdf"
# For doc coordinates, (0, 0) is the bottom left corner


# TODO: want a fn that allows you to draw a multi-page table that has repeated headers
# so the first page the table is introduced has a different header e.g. FOO
# and all others are like FOO (continued)
# Reportlab doesn't support this well, but you can hack this by creating separate
# pdf files for each table. Then merge everything together and the header and footer
# This will be a generic fn to reuse will all the tables in the report
# def table_with_headers() -> :


def default_shapes(canvas, doc: BaseDocTemplate):
    canvas.saveState()
    bottom_line_h = inch * 1
    top_line_h = PHEIGHT - inch * 1
    canvas.line(0, bottom_line_h, PWIDTH, bottom_line_h)
    canvas.line(0, top_line_h, PWIDTH, top_line_h)
    canvas.restoreState()


# Using Paragraphs to format text for word wrap
column_style: ParagraphStyle = ParagraphStyle("cols", fontSize=9)
numeric_style: ParagraphStyle = ParagraphStyle("nums", fontSize=11)
text_style: ParagraphStyle = ParagraphStyle("data", fontSize=10)
data_styles: dict = {1: numeric_style, 2: numeric_style, None: text_style}
ncols: int = len(res.columns)

# # Dummy styles for now
col_styles = rp.style_cells(
    (0, 0),
    ncols,
    1,
    textcolor=colors.red,
    fontsize=13,
    underline=(3, colors.black),
    background=colors.lightgrey,
)
style_all = rp.style_cells((0, 1), background=colors.lightcyan, valign="TOP")
test_table = rp.ReportElement(
    "/home/shannc/Bio_SDD/chula-stem/test2.pdf",
    {"header_pos": (inch * 5, A4[1] - cm)},
    header_first="Relevant Variants",
    header_later="Relevant Variants (continued)",
)
test_table.add_table(
    res,
    data_styles,
    column_style,
    style_all,
    col_styles,
    col_widths=list(VTABLE_COL_WIDTHS.values()),
)
test_table.add_decorator(default_shapes)
test_table.build()
