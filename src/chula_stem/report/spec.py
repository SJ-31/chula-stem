from dataclasses import dataclass
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from chula_stem.report.utils import style_cells
from reportlab.lib import colors

STYLES = getSampleStyleSheet()

FONT = "Helvetica"
BOLD_FONT = "Helvetica-Bold"
NUMERIC_STYLE: ParagraphStyle = ParagraphStyle("nums", fontSize=11)
TEXT_STYLE: ParagraphStyle = ParagraphStyle("text", fontSize=10)
TEXT_STYLE_SMALL: ParagraphStyle = ParagraphStyle("text_small", fontSize=8)

# Use paragraph styles to specify font parameters (e.g. color, size, font face)
# and not the cell styles


@dataclass
class URL:
    dbvar = "https://www.ncbi.nlm.nih.gov/dbvar/variants"
    pubchem = "https://pubchem.ncbi.nlm.nih.gov/compound"


@dataclass
class Widths:
    sv_snp = {
        "Gene": 80,
        "Variant Allele Frequency": 30,
        "Variant Read Support": 40,
        "Locus": 80,
        "HGVS": 80,
        "Database Name": 90,
        "Variant Type": 90,
        "ClinVar": 90,
    }
    cnv = {
        "Locus": 160,
        "CNV Type": 50,
        "Estimated Copy Number": 70,
        "Known/predicted Dosage-sensitive Genes": 90,
        "All Genes": 90,
        "ClinGen": 50,
        "Database/Study Records": 70,
    }
    reference_table = {"id": 30, "text": 500}
    repeat = {
        "Locus": 130,
        "Repeat Unit": 80,
        "Repeat Number": 50,
        "Affected Gene": 90,
        "ClinGen": 50,
        "Database/Study Records": 100,
    }
    signature = {
        "Signature": 60,
        "Count": 60,
        "Frequency": 60,
        "Collection": 60,
        "Proposed Aetiology": 150,
    }
    therapy = {
        "Therapy": 100,
        "PubChemId": 70,
        "Evidence category": 90,
        "Evidence in sample": 90,
        "Relevant cancers": 100,
        "Study": 60,
        "Database source": 70,
    }


@dataclass
class Rename:
    snp = {
        "SYMBOL": "Gene",
        "VAF": "VAF (%)",
        "Alt_depth": "Variant Read Support",
        "Loc": "Locus",
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
        "Known or predicted dosage-sensitive genes": "Known/predicted Dosage-sensitive Genes",
        "All protein coding genes": "All Genes",
        "ClinGen_report": "ClinGen",
        "source": "Database/Study Records",
    }
    therapy = {
        "therapies": "Therapy",
        "PubChemId": "PubChemId",
        "type": "Evidence category",
        "gene": "Evidence in sample",
        "disease": "Relevant cancers",
        "source": "Study source(s)",
        "db_link": "Database source(s)",
    }
    repeat = {
        "repeat_unit_bases": "Repeat Unit",
        "repeat_times": "Repeat Number",
        "gene_name": "Affected Gene",
        "ClinGen_report": "ClinGen",
        "source": "Database/Study Records",
    }


def snp_style():
    cell_pstyles: dict = {
        None: TEXT_STYLE_SMALL,
    }
    header_pstyles: dict = {
        None: ParagraphStyle("header", fontSize=9, fontName=BOLD_FONT),
        2: ParagraphStyle("header_small", fontSize=8, fontName=BOLD_FONT),
        1: ParagraphStyle("header_small", fontSize=8, fontName=BOLD_FONT),
    }
    cell_styles = style_cells((0, 1), background=colors.whitesmoke, valign="TOP")
    header_styles = style_cells(
        (0, 0),
        8,
        1,
        underline=(3, colors.black),
        background=colors.lightgrey,
    )
    return {
        "cell_pstyles": cell_pstyles,
        "header_pstyles": header_pstyles,
        "cell_styles": cell_styles,
        "header_styles": header_styles,
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
    cell_pstyles: dict = {None: TEXT_STYLE_SMALL}
    cell_styles = style_cells((0, 1), background=colors.whitesmoke, valign="TOP")
    header_styles = style_cells(
        (0, 0),
        7,
        1,
        underline=(3, colors.black),
        background=colors.lightgrey,
    )
    return {
        "cell_pstyles": cell_pstyles,
        "header_pstyles": ParagraphStyle("cols", fontSize=9, fontName=BOLD_FONT),
        "cell_styles": cell_styles,
        "header_styles": header_styles,
        "col_widths": list(Widths.therapy.values()),
    }


def cnv_style():
    cell_pstyles: dict = {0: NUMERIC_STYLE, 2: NUMERIC_STYLE, None: TEXT_STYLE_SMALL}
    cell_styles = style_cells(
        (0, 1), fontname=FONT, background=colors.whitesmoke, valign="TOP"
    )
    header_styles = style_cells(
        (0, 0),
        7,
        1,
        textcolor=colors.red,
        underline=(3, colors.black),
        background=colors.lightgrey,
    )
    return {
        "cell_pstyles": cell_pstyles,
        "header_pstyles": ParagraphStyle("cols", fontSize=9, fontName=BOLD_FONT),
        "cell_styles": cell_styles,
        "header_styles": header_styles,
        "col_widths": list(Widths.cnv.values()),
    }


def signature_style():
    cell_pstyles: dict = {1: NUMERIC_STYLE, None: TEXT_STYLE}
    cell_styles = style_cells(
        (0, 1), fontname=FONT, background=colors.whitesmoke, valign="TOP"
    )
    header_styles = style_cells(
        (0, 0),
        6,
        1,
        textcolor=colors.red,
        fontname=BOLD_FONT,
        underline=(3, colors.black),
        background=colors.lightgrey,
    )
    return {
        "cell_pstyles": cell_pstyles,
        "header_pstyles": ParagraphStyle("cols", fontSize=9),
        "cell_styles": cell_styles,
        "header_styles": header_styles,
        "col_widths": list(Widths.signature.values()),
    }


def repeat_style():
    cell_pstyles: dict = {2: NUMERIC_STYLE, None: TEXT_STYLE}
    cell_styles = style_cells((0, 1), background=colors.whitesmoke, valign="TOP")
    header_styles = style_cells(
        (0, 0),
        6,
        1,
        textcolor=colors.red,
        underline=(3, colors.black),
        background=colors.lightgrey,
    )
    return {
        "cell_pstyles": cell_pstyles,
        "header_pstyles": ParagraphStyle("cols", fontSize=9, fontName=BOLD_FONT),
        "cell_styles": cell_styles,
        "header_styles": header_styles,
        "col_widths": list(Widths.repeat.values()),
    }
