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
    "table_title_style": ParagraphStyle(
        "headline", fontSize=14, fontName=BOLD_FONT, textColor=HexColor("#549f9f")
    ),
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


@dataclass
class URL:
    dbvar = "https://www.ncbi.nlm.nih.gov/dbvar/variants"
    pubchem = "https://pubchem.ncbi.nlm.nih.gov/compound"


@dataclass
class Widths:
    available_width = A4[0] - 4 * cm
    sv_snp = {
        "Locus": available_width * 0.13,
        "Variant Allele Frequency": available_width * 0.08,
        "Variant Read Support": available_width * 0.07,
        "Gene": available_width * 0.12,
        "HGVS": available_width * 0.14,
        "Database Name": available_width * 0.16,
        "Variant Type": available_width * 0.2,
        "ClinVar": available_width * 0.1,
    }
    cnv = {
        "Locus": available_width * 0.21,
        "CNV Type": available_width * 0.05,
        "Estimated Copy Number": available_width * 0.1,
        "Known/predicted Dosage-sensitive Genes": available_width * 0.15,
        "All Genes": available_width * 0.29,
        "ClinGen": available_width * 0.08,
        "Database/Study Records": available_width * 0.12,
    }
    reference_table = {"id": 30, "text": 500}
    repeat = {
        "Locus": available_width * 0.21,
        "Repeat Unit": available_width * 0.12,
        "Repeat Number": available_width * 0.08,
        "Affected Gene": available_width * 0.12,
        "ClinGen": available_width * 0.08,
        "Database/Study Records": available_width * 0.39,
    }
    signature = {
        "Signature": available_width * 0.12,
        "Count": available_width * 0.15,
        "Frequency": available_width * 0.1,
        "Collection": available_width * 0.1,
        "Proposed Aetiology": available_width * 0.53,
    }
    therapy = {
        "Therapy": available_width * 0.19,
        "PubChemId": available_width * 0.1,
        "Evidence category": available_width * 0.16,
        "Evidence in sample": available_width * 0.15,
        "Relevant cancers": available_width * 0.2,
        "Study": available_width * 0.1,
        "Database source": available_width * 0.1,
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
        "VAF": "VAF (%)",
        "Loc": "Locus",
        "HGVSc": "HGVS",
        "Existing_variation": "Database Name",
        "Consequence": "Variant Type",
        "SVTYPE": "SV Class",
        "CLIN_SIG": "ClinVar",
    }
    cnv = {
        "Type": "CNV Type",
        "Known or predicted dosage-sensitive genes": "Dosage-sensitive Genes",
        "All protein coding genes": "All Genes",
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
        "cell_styles": STYLE["cell_style"](Widths.sv_snp),
        "header_styles": STYLE["header_style"](Widths.sv_snp),
        "col_widths": list(Widths.sv_snp.values()),
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
