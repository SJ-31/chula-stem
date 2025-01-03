import os
from collections import Counter
from pathlib import Path
from typing import Callable

import polars as pl
import polars.selectors as cs
import pytest
from chula_stem.callset_qc import IMPACT_MAP
from chula_stem.report.format import msisensor_pro_fmt

small_path = "/home/shannc/Bio_SDD/chula-stem/tests/vep/7-patient_10-VEP_small_1.tsv"
sv_pat = "/home/shannc/Bio_SDD/chula-stem/tests/vep/sv2.tsv"

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


@pytest.mark.skip(reason="Done")
def test_cnvkit_pdf():
    from chula_stem.report import VariantCallingReport

    cnr = "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.cnr"
    cns = "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.call.cns"
    VariantCallingReport.plot_cnvkit(cnr, cns)


# @pytest.mark.skip(reason="Done")
def test_full():
    from chula_stem.report.variant_calling_report import VariantCallingReport

    meta = {
        "patient_name": "John Doe",
        "gender": "Male",
        "patient_birthdate": "1990-05-14",
        "physician": "Dr. Alice Smith",
        "diagnosis": "Hypertension",
        "hospital": "Green Valley Medical Center",
        "id": "SMP-123456",
        "collection_date": "2024-12-01",
        "sample_type": "Blood",
    }
    R = VariantCallingReport(
        "/home/shannc/Bio_SDD/chula-stem/tests/report_full.pdf",
        metadata=meta,
        civic_cache="/home/shannc/Bio_SDD/chula-stem/tests/civic.json",
        pandrugs2_cache="/home/shannc/Bio_SDD/chula-stem/tests/pandrugs2.json",
        vep_small=small_path,
        vep_sv=sv_pat,
        classify_cnv="/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-classify_cnv-CR.tsv",
        facets="/home/shannc/Bio_SDD/chula-stem/tests/5-sample2-Facets/5-patient_10-Facets_hisens.rds",
        cnvkit_cns="/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.call.cns",
        cnvkit_cnr="/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.cnr",
        msisensor_pro="/home/shannc/Bio_SDD/chula-stem/tests/msisensor/4-null-CR.tsv",
        sigprofiler="/home/shannc/Bio_SDD/chula-stem/tests/6-null-SigProfilerAssignment/Activities/Assignment_Solution_Activities.txt",
        cosmic_reference="/home/shannc/Bio_SDD/chula-stem/nextflow/config/cosmic_signatures_v3.4-2024-12-26.csv",
        tmpdir="/home/shannc/Bio_SDD/chula-stem/tests/report_tmp",
        plot=False,
    )
    R.build()
    R.merge()


@pytest.mark.skip(reason="Done")
def test_format_classify():
    from chula_stem.report.format import classify_cnv_fmt

    classify_cnv = (
        "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-classify_cnv-CR.tsv"
    )
    facets = "/home/shannc/Bio_SDD/chula-stem/tests/5-sample2-Facets/5-patient_10-Facets_hisens.rds"
    cnvkit = "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.call.cns"
    classify_cnv_fmt(classify_cnv, facets, cnvkit)


@pytest.mark.skip(reason="Done")
def test_format_vep_sv():
    from chula_stem.report.format import vep_fmt

    try:
        os.mkdir("/home/shannc/Bio_SDD/chula-stem/tests/vep_format_sv")
    except:
        pass
    vep_fmt(sv_pat, variant_class="sv")
    vep_fmt(small_path, variant_class="small")


spec: list = [
    {
        "file": "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-classify_cnv-CR.tsv",
        "is_list": True,
        "is_list_separator": ",",
        "column": "All protein coding genes",
    },
    {
        "file": "/home/shannc/Bio_SDD/chula-stem/tests/msisensor/5-P1-Msisensor_unstable.tsv",
        "column": "gene_name",
    },
    {
        "file": "/home/shannc/Bio_SDD/chula-stem/tests/vep/sv.tsv",
        "column": "SYMBOL",
    },
]
