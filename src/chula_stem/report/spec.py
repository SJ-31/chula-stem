from dataclasses import dataclass

from reportlab.lib import colors
from reportlab.lib.colors import HexColor
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm, mm
from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus import BaseDocTemplate

from chula_stem.report.utils import alternating_bg, style_cells

FONT = "Helvetica"
BOLD_FONT = "Helvetica-Bold"
NUMERIC_STYLE: ParagraphStyle = ParagraphStyle("nums", fontSize=11)
TEXT_STYLE: ParagraphStyle = ParagraphStyle("text", fontSize=10)
TEXT_STYLE_SMALL: ParagraphStyle = ParagraphStyle("text_small", fontSize=8)

# * Variant calling report style 1
# <2025-01-02 Thu>
# Author: Manoon Chongwattananukul
# Usd first for whole-exome sequencing pipeline
STYLE: dict = {
    "header_bg": "#559f9f",
    "cell1": "#c3e7eb",
    "cell2": "#87cfd6",
    "cell_font_size": 8,
    "cell_font_size_small": 6,
    "header_pstyle": ParagraphStyle(
        "header", fontSize=6.5, fontName=BOLD_FONT, textColor=colors.white, alignment=1
    ),
    "toc_title_style": ParagraphStyle("toc", fontSize=14),
    "table_title_style": ParagraphStyle(
        "headline", fontSize=14, fontName=BOLD_FONT, textColor=HexColor("#549f9f")
    ),
    "toc_text": lambda x: ParagraphStyle(
        "text", fontName=FONT, fontSize=8, firstLineIndent=(x - 1) * 20
    ),
    "toc_num": ParagraphStyle("num", fontName=BOLD_FONT, fontSize=8, alignment=2),
}
STYLE["cell_pstyle"] = ParagraphStyle("text", fontSize=6.5)
STYLE["cell_style"] = lambda x: alternating_bg(
    len(x),
    STYLE["cell1"],
    STYLE["cell2"],
    offset=1,
    valign="TOP",
    grid=(1, colors.white),
)
STYLE["header_style"] = lambda x: style_cells(
    (0, 0),
    len(x),
    1,
    underline=(1, colors.white),
    background=STYLE["header_bg"],
    valign="TOP",
    align="CENTER",
    grid=(1, colors.white),
)
STYLE["toc_cell_bg"] = "#c3e7ea"
STYLE["title_style"] = ParagraphStyle(
    "tstyle",
    textColor=colors.white,
    fontSize=25,
    fontName=FONT,
    leading=25,
    alignment=1,
)
STYLE["detail_style"] = ParagraphStyle(
    "dstyle", fontSize=8, fontName=FONT, firstLineIndent=10
)


def decorator(c: Canvas, d: BaseDocTemplate) -> None:
    c.saveState()
    x = 2 * cm
    xend = A4[0] - 2 * cm
    y = A4[1] - 1.9 * cm
    c.setLineWidth(0.9)
    c.setStrokeColor(HexColor("#86cfd5"))
    c.line(x, y, xend, y)
    c.restoreState()


STYLE["table_decorator"] = decorator
AVAILABLE_WIDTH = A4[0] - 4 * cm


@dataclass
class URL:
    dbvar = "https://www.ncbi.nlm.nih.gov/dbvar/variants"
    pubchem = "https://pubchem.ncbi.nlm.nih.gov/compound"


@dataclass
class Widths:
    snp = {
        "Locus": AVAILABLE_WIDTH * 0.13,
        "Variant Allele Frequency": AVAILABLE_WIDTH * 0.08,
        "Variant Read Support": AVAILABLE_WIDTH * 0.07,
        "Gene": AVAILABLE_WIDTH * 0.12,
        "HGVS": AVAILABLE_WIDTH * 0.14,
        "Database Name": AVAILABLE_WIDTH * 0.16,
        "Variant Type": AVAILABLE_WIDTH * 0.2,
        "ClinVar": AVAILABLE_WIDTH * 0.1,
    }
    sv = {
        "Locus": AVAILABLE_WIDTH * 0.13,
        "Gene": AVAILABLE_WIDTH * 0.12,
        "HGVS": AVAILABLE_WIDTH * 0.14,
        "Database Name": AVAILABLE_WIDTH * 0.16,
        "Variant Type": AVAILABLE_WIDTH * 0.2,
        "SV Class": AVAILABLE_WIDTH * 0.15,
        "ClinVar": AVAILABLE_WIDTH * 0.1,
    }
    cnv = {
        "Locus": AVAILABLE_WIDTH * 0.21,
        "CNV Type": AVAILABLE_WIDTH * 0.05,
        "Estimated Copy Number": AVAILABLE_WIDTH * 0.1,
        "ClinVar": AVAILABLE_WIDTH * 0.10,
        "All Genes": AVAILABLE_WIDTH * 0.34,
        "ClinGen": AVAILABLE_WIDTH * 0.08,
        "Database/Study Records": AVAILABLE_WIDTH * 0.12,
    }
    reference_table = {"id": 30, "text": 500}
    repeat = {
        "Locus": AVAILABLE_WIDTH * 0.21,
        "Repeat Unit": AVAILABLE_WIDTH * 0.12,
        "Repeat Number": AVAILABLE_WIDTH * 0.08,
        "Affected Gene": AVAILABLE_WIDTH * 0.12,
        "ClinGen": AVAILABLE_WIDTH * 0.08,
        "Database/Study Records": AVAILABLE_WIDTH * 0.39,
    }
    signature = {
        "Signature": AVAILABLE_WIDTH * 0.12,
        "Count": AVAILABLE_WIDTH * 0.18,
        "Frequency": AVAILABLE_WIDTH * 0.1,
        "Collection": AVAILABLE_WIDTH * 0.1,
        "Proposed Aetiology": AVAILABLE_WIDTH * 0.50,
    }
    therapy = {
        "Therapy": AVAILABLE_WIDTH * 0.19,
        "PubChemId": AVAILABLE_WIDTH * 0.1,
        "Evidence category": AVAILABLE_WIDTH * 0.16,
        "Evidence in sample": AVAILABLE_WIDTH * 0.15,
        "Relevant cancers": AVAILABLE_WIDTH * 0.2,
        "Study": AVAILABLE_WIDTH * 0.1,
        "Database source": AVAILABLE_WIDTH * 0.1,
    }


@dataclass
class Rename:
    snp = {
        "Loc": "Locus",
        "VAF": "VAF (%)",
        "Alt_depth": "VRS",
        "SYMBOL": "Gene",
        "HGVSc": "HGVS",
        "Existing_variation": "Database Name",
        "Consequence": "Variant Type",
        "CLIN_SIG": "ClinVar",
    }
    sv = {
        "SYMBOL": "Gene",
        "Loc": "Locus",
        "HGVSc": "HGVS",
        "Existing_variation": "Database Name",
        "Consequence": "Variant Type",
        "SVTYPE": "SV Class",
        "CLIN_SIG": "ClinVar",
    }
    cnv = {
        "Type": "CNV Type",
        "Classification": "ClinVar",
        "All protein coding genes": "Genes in Region",
        "ClinGen_report": "ClinGen",
        "source": "Source",
    }
    therapy = {
        "therapies": "Therapy",
        "PubChemId": "PubChemId",
        "type": "Evidence category",
        "gene": "Evidence in sample",
        "disease": "Relevant cancers",
        "source": "Study source",
        "db_link": "Database source",
    }
    repeat = {
        "repeat_unit_bases": "Repeat Unit",
        "repeat_times": "Repeat Number",
        "gene_name": "Affected Gene",
        "ClinGen_report": "ClinGen",
        "source": "Source",
    }


def snp_style():
    cell_pstyles = {
        None: STYLE["cell_pstyle"],
        1: ParagraphStyle("text", fontSize=6.5, alignment=1),
        2: ParagraphStyle("text", fontSize=6.5, alignment=1),
    }
    return {
        "cell_pstyles": cell_pstyles,
        "header_pstyles": STYLE["header_pstyle"],
        "cell_styles": STYLE["cell_style"](Widths.snp),
        "header_styles": STYLE["header_style"](Widths.snp),
        "col_widths": list(Widths.snp.values()),
    }


def sv_style():
    cell_pstyles = {None: STYLE["cell_pstyle"]}
    return {
        "cell_pstyles": cell_pstyles,
        "header_pstyles": STYLE["header_pstyle"],
        "cell_styles": STYLE["cell_style"](Widths.sv),
        "header_styles": STYLE["header_style"](Widths.sv),
        "col_widths": list(Widths.sv.values()),
    }


def reference_list_style():
    cell_pstyles: dict = {
        None: TEXT_STYLE,
    }
    cell_styles = style_cells((0, 1), background=colors.white, valign="TOP")
    return {
        "cell_pstyles": cell_pstyles,
        "header_pstyles": None,
        "cell_styles": cell_styles,
        "header_styles": None,
        "col_widths": list(Widths.reference_table.values()),
        "repeat_header": False,
    }


def therapy_style():
    cell_pstyles = {
        None: STYLE["cell_pstyle"],
        1: ParagraphStyle("text", fontSize=6.5, alignment=1),
    }
    return {
        "cell_pstyles": cell_pstyles,
        "header_pstyles": STYLE["header_pstyle"],
        "cell_styles": STYLE["cell_style"](Widths.therapy),
        "header_styles": STYLE["header_style"](Widths.therapy),
        "col_widths": list(Widths.therapy.values()),
    }


def cnv_style():
    cell_pstyles = {
        None: STYLE["cell_pstyle"],
        1: ParagraphStyle("text", fontSize=6.5, alignment=1),
        2: ParagraphStyle("text", fontSize=6.5, alignment=1),
    }
    return {
        "cell_pstyles": cell_pstyles,
        "header_pstyles": STYLE["header_pstyle"],
        "cell_styles": STYLE["cell_style"](Widths.cnv),
        "header_styles": STYLE["header_style"](Widths.cnv),
        "col_widths": list(Widths.cnv.values()),
    }


def signature_style():
    cell_pstyles = {
        None: STYLE["cell_pstyle"],
        1: ParagraphStyle("text", fontSize=6.5, alignment=1),
        2: ParagraphStyle("text", fontSize=6.5, alignment=1),
        3: ParagraphStyle("text", fontSize=6.5, alignment=1),
    }
    return {
        "cell_pstyles": cell_pstyles,
        "header_pstyles": STYLE["header_pstyle"],
        "cell_styles": STYLE["cell_style"](Widths.signature),
        "header_styles": STYLE["header_style"](Widths.signature),
        "col_widths": list(Widths.signature.values()),
    }


def repeat_style():
    cell_pstyles = {
        None: STYLE["cell_pstyle"],
        2: ParagraphStyle("text", fontSize=6.5, alignment=1),
    }
    return {
        "cell_pstyles": cell_pstyles,
        "header_pstyles": STYLE["header_pstyle"],
        "cell_styles": STYLE["cell_style"](Widths.repeat),
        "header_styles": STYLE["header_style"](Widths.repeat),
        "col_widths": list(Widths.repeat.values()),
    }
