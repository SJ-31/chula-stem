import os
from collections import Counter
from pathlib import Path
from typing import Callable

import polars as pl
import polars.selectors as cs
from chula_stem.callset_qc import IMPACT_MAP

small_path = "/home/shannc/Bio_SDD/chula-stem/tests/vep/7-patient_10-VEP_small_1.tsv"
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


def default_shapes(canvas, doc: BaseDocTemplate):
    canvas.saveState()
    bottom_line_h = inch * 1
    top_line_h = PHEIGHT - inch * 1
    canvas.line(0, bottom_line_h, PWIDTH, bottom_line_h)
    canvas.line(0, top_line_h, PWIDTH, top_line_h)
    canvas.restoreState()


def test_full():
    from chula_stem.report import ResultsReport

    R = ResultsReport(
        "/home/shannc/Bio_SDD/chula-stem/report_full.pdf",
        civic_cache="/home/shannc/Bio_SDD/chula-stem/tests/civic.json",
        pandrugs2_cache="/home/shannc/Bio_SDD/chula-stem/tests/pandrugs2.json",
        vep_small=small_path,
        vep_sv=sv_pat,
        classify_cnv="/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-classify_cnv-CR.tsv",
        facets="/home/shannc/Bio_SDD/chula-stem/tests/5-sample2-Facets/5-patient_10-Facets_hisens.rds",
        cnvkit="/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.call.cns",
        tmpdir="/home/shannc/Bio_SDD/chula-stem/tests/report_tmp",
    )
    R.build()


def test_format_classify():
    from chula_stem.report.format import classify_cnv_fmt

    classify_cnv = (
        "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-classify_cnv-CR.tsv"
    )
    facets = "/home/shannc/Bio_SDD/chula-stem/tests/5-sample2-Facets/5-patient_10-Facets_hisens.rds"
    cnvkit = "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.call.cns"
    classify_cnv_fmt(classify_cnv, facets, cnvkit)
