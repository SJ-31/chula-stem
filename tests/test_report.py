from pathlib import Path

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
import polars.selectors as cs
from reportlab.lib import colors, units
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.pdfgen import canvas
from reportlab.platypus import LongTable, Paragraph, SimpleDocTemplate, Spacer

# (0, 0) is the bottom left of the page by default

# Flowables
# - Flowable.drawOn(canvas, x, y)

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

# Using Paragraphs to format text for word wrap
column_style: ParagraphStyle = ParagraphStyle("cols", fontSize=9)
numeric_style: ParagraphStyle = ParagraphStyle("nums", fontSize=11)
text_style: ParagraphStyle = ParagraphStyle("data", fontSize=10)
data_styles: dict = {1: numeric_style, 2: numeric_style, None: text_style}

to_data = rp.add_pstyles(res, data_styles)
cols = rp.add_pstyles(res.columns, column_style)
to_data.insert(0, cols)

# repeatRows=1 repeats the first row at every split
T = LongTable(to_data, colWidths=list(VTABLE_COL_WIDTHS.values()), repeatRows=1)
ncols: int = len(res.columns)

# Dummy styles for now
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

# Use a doc template to enable the table to be split across multiple pages
pdf = SimpleDocTemplate(
    "/home/shannc/Bio_SDD/chula-stem/test.pdf", rightMargin=3, leftMargin=3
)
T.setStyle(col_styles + style_all)
loc: tuple = (10, 10)
# splits = T.splitOn(pdf, 500, 100)
elements = []
elements.append(T)


def _header_footer(canvas, doc):
    # Save the state of our canvas so we can draw on it
    canvas.saveState()
    styles = getSampleStyleSheet()

    # Header
    header = Paragraph(
        "This is a multi-line header.  It goes on every page.   " * 5,
        styles["Normal"],
    )
    w, h = header.wrap(doc.width, doc.topMargin)
    header.drawOn(canvas, doc.leftMargin, doc.height + doc.topMargin - h)

    # Footer
    footer = Paragraph(
        "This is a multi-line footer.  It goes on every page.   " * 5,
        styles["Normal"],
    )
    w, h = footer.wrap(doc.width, doc.bottomMargin)
    footer.drawOn(canvas, doc.leftMargin, h)

    # Release the canvas
    canvas.restoreState()


pdf.build(elements, onFirstPage=_header_footer, onLaterPages=_header_footer)
