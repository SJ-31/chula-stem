from dataclasses import dataclass
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from chula_stem.report.utils import style_cells
from reportlab.lib import colors

STYLES = getSampleStyleSheet()


@dataclass
class URL:
    dbvar = "https://www.ncbi.nlm.nih.gov/dbvar/variants"
    pubchem = "https://pubchem.ncbi.nlm.nih.gov/compound"


@dataclass
class Widths:
    sv_snp = {
        "Gene": 60,
        "Variant Allele Frequency": 55,
        "Variant Read Support": 45,
        "Locus": 80,
        "HGVS": 90,
        "Database Name": 90,
        "Variant Type": 90,
        "ClinVar": 100,
    }
    cnv = {
        "Locus": 100,
        "CNV Type": 50,
        "Estimated Copy Number": 40,
        "Known/predicted Dosage-sensitive Genes": 90,
        "All Genes": 100,
        "ClinGen": 50,
        "Database/Study Records": 50,
    }
    repeat = {
        "Locus": 50,
        "Repeat Unit": 50,
        "Repeat Number": 50,
        "Affected Gene": 50,
        "ClinGen": 50,
        "Database/Study Records": 50,
    }
    therapy = {
        "Therapy": 90,
        "PubChemId": 90,
        "Evidence category": 70,
        "Evidence in sample": 90,
        "Relevant cancers": 100,
        "Study": 100,
        "Database source": 90,
    }


@dataclass
class Rename:
    sv_snp = {
        "SYMBOL": "Gene",
        "VAF": "Variant Allele Frequency",
        "Alt_depth": "Variant Read Support",
        "Loc": "Locus",
        "HGVSc": "HGVS",
        "Existing_variation": "Database Name",
        "Consequence": "Variant Type",
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
        "Gene": "Evidence in sample",
        "disease": "Relevant cancers",
        "source": "Study source",
        "db_link": "Database source",
    }
    repeat = {
        "repeat_unit_bases": "Repeat Unit",
        "repeat_times": "Repeat Number",
        "gene_name": "Affected Gene",
        "ClinGen_report": "ClinGen",
        "source": "Database/Study Records",
    }


NUMERIC_STYLE: ParagraphStyle = ParagraphStyle("nums", fontSize=11)
TEXT_STYLE: ParagraphStyle = ParagraphStyle("data", fontSize=10)


def sv_snp_style():
    cell_pstyles: dict = {
        1: NUMERIC_STYLE,
        2: NUMERIC_STYLE,
        None: TEXT_STYLE,
    }
    cell_styles = style_cells((0, 1), background=colors.lightcyan, valign="TOP")
    header_styles = style_cells(
        (0, 0),
        8,
        1,
        textcolor=colors.red,
        underline=(3, colors.black),
        background=colors.lightgrey,
    )
    return {
        "cell_pstyles": cell_pstyles,
        "header_pstyles": ParagraphStyle("cols", fontSize=9),
        "cell_styles": cell_styles,
        "header_styles": header_styles,
        "col_widths": list(Widths.sv_snp.values()),
    }


def therapy_style():
    cell_pstyles: dict = {
        2: NUMERIC_STYLE,
        None: TEXT_STYLE,
    }
    cell_styles = style_cells((0, 1), background=colors.lightcyan, valign="TOP")
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
        "header_pstyles": ParagraphStyle("cols", fontSize=9),
        "cell_styles": cell_styles,
        "header_styles": header_styles,
        "col_widths": list(Widths.therapy.values()),
    }


def cnv_style():
    cell_pstyles: dict = {0: NUMERIC_STYLE, 2: NUMERIC_STYLE, None: TEXT_STYLE}
    cell_styles = style_cells((0, 1), background=colors.lightcyan, valign="TOP")
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
        "header_pstyles": ParagraphStyle("cols", fontSize=9),
        "cell_styles": cell_styles,
        "header_styles": header_styles,
        "col_widths": list(Widths.cnv.values()),
    }


def repeat_style():
    cell_pstyles: dict = {2: NUMERIC_STYLE, None: TEXT_STYLE}
    cell_styles = style_cells((0, 1), background=colors.lightcyan, valign="TOP")
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
        "header_pstyles": ParagraphStyle("cols", fontSize=9),
        "cell_styles": cell_styles,
        "header_styles": header_styles,
        "col_widths": list(Widths.repeat.values()),
    }
