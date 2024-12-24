import os
from typing import Callable

import polars as pl
import polars.selectors as cs
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
)
from chula_stem.report.spec import (
    Rename,
    Widths,
)

STYLES = getSampleStyleSheet()

import chula_stem.report.format as fr


def filter_format_vep(input: str, sep="\t"):
    wanted_cols = [
        "Gene",
        "VAF",
        "Alt_depth",
        "Loc",
        "SYMBOL",
        "Consequence",
        "CLIN_SIG",
        "Known Variants",
        "HGVSg",
    ]
    to_rename: dict = {
        "Alt_depth": "Variant Allele Depth",
        "VAF": "Variant Allele Fraction",
        "Loc": "Location",
        "CLIN_SIG": "ClinVar",
        "HGVSg": "HGVS",
        "Existing_variation": "Known Variants",
        "Consequence": "Variant Effect",
    }
    df = (
        pl.read_csv(input, separator=sep)
        .with_columns(
            pl.col("Existing_variation").map_elements(
                lambda x: x.replace("&", ", "), return_dtype=pl.String
            )
        )
        .select(wanted_cols)
        .rename(to_rename)
    )
    return df


## * Report formatter


def draw_paragraph(
    text: str,
    pos: tuple,
    style: ParagraphStyle,
    canvas: Canvas,
    doc: BaseDocTemplate,
):
    canvas.saveState()
    para = Paragraph(text, style)
    para.wrap(doc.width, doc.topMargin)
    para.drawOn(canvas, pos[0], pos[1])
    canvas.restoreState()


def style_cells(
    start: tuple,
    ncols: int = 0,
    nrows: int = 0,
    fontname: str = "",
    fontsize: str = "",
    align: str = "",
    background: str = "",
    valign: str = "",
    textcolor=None,
    underline: tuple = (),
) -> list:
    """
    Coordinates for table style are given as (column, row)
    """
    parameter_map: dict = {
        fontname: "FONTNAME",
        fontsize: "FONTSIZE",
        textcolor: "TEXTCOLOR",
        background: "BACKGROUND",
        align: "ALIGN",
        valign: "VALIGN",
        underline: "LINEBELOW",
    }
    if ncols and nrows:
        end: tuple = start[0] + ncols - 1, start[1] + nrows - 1
    else:
        end = (-1, -1)

    def style_helper(format: str, value) -> tuple:
        if isinstance(value, tuple):
            return (format, start, end) + value
        else:
            return (format, start, end, value)

    styles: list = []
    for param, name in parameter_map.items():
        if param:
            styles.append(style_helper(name, param))
    return styles


def add_pstyles(
    data: pl.DataFrame | list, style: ParagraphStyle | dict[int, ParagraphStyle]
) -> list:
    """Helper for adding paragraph styles to lists or data in dfs

    :param style: A single style which is then applied to all data.
    Alternatively, a map of column_index -> style specifying styles to
    apply to specific columns. A key for 'None' is the default and applied to columns
    not explicitly given

    :returns:
    """
    if isinstance(style, dict):
        if not style.get(None):
            raise ValueError("A default key `None` must be provided!")

        style_fn = lambda x, index=None: Paragraph(str(x), style.get(index))
    else:
        style_fn = lambda x, index=None: Paragraph(str(x), style)

    if isinstance(data, pl.DataFrame):
        return [
            [Paragraph(str(s), style.get(i)) for i, s in enumerate(row)]
            for row in data.iter_rows()
        ]
    else:
        return [style_fn(row, i) for i, row in enumerate(data)]


## ** Report class


class ResultsReport:
    def __init__(
        self,
        filename: str,
        pandrugs2_cache: str = "",
        civic_cache: str = "",
        vep_small: str = "",
        vep_sv: str = "",
        classify_cnv: str = "",
        facets: str = "",
        cnvkit: str = "",
        tmpdir: str = "temp",
    ) -> None:
        self.civic_cache = civic_cache
        self.pandrugs2_cache = pandrugs2_cache
        self.data: dict = {"relevant": {}, "nonrelevant": {}, "all": {}}
        self.tmpdir = tmpdir

        try:
            os.makedirs(self.tmpdir)
        except FileExistsError:
            print("WARNING: directory exists")
        os.chdir(self.tmpdir)
        calls = [
            lambda: fr.vep_fmt(vep_small, tmpdir, "small", civic_cache, pandrugs2_cache),
            lambda: fr.vep_fmt(vep_sv, tmpdir, "sv", civic_cache, pandrugs2_cache),
            lambda: fr.classify_cnv_fmt(classify_cnv, facets, cnvkit),
        ]
        for type, fn_call in zip(
            ["small", "sv", "cnv"],
            calls,
        ):
            all, relevant, nonrelevant = fn_call()
            self.data["relevant"][type] = relevant
            self.data["nonrelevant"][type] = nonrelevant
            self.data["all"][type] = all

    @staticmethod
    def build_table(
        spec: dict,
        table_styles: dict,
        data: pl.DataFrame,
        first: str,
        later: str,
        decorator: Callable | None,
        filename: str,
    ) -> None:
        spec = {**spec, "header_first": first, "header_later": later}
        R = ReportElement(filename, spec)
        R.add_table(data, **table_styles)
        R.add_decorator(decorator)
        R.build()

    def build(self) -> None:
        """Create pdf files for all report elements individually, concatenate them
        and add header (time + page number)
        """
        table_spec = {"header_pos": (A4[0] - 5 * inch, A4[1] - inch)}
        numeric_style: ParagraphStyle = ParagraphStyle("nums", fontSize=11)
        text_style: ParagraphStyle = ParagraphStyle("data", fontSize=10)

        cell_pstyles: dict = {
            1: numeric_style,
            2: numeric_style,
            None: text_style,
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
        # Variant table styles

        vtable_styles: dict = {
            "cell_pstyles": cell_pstyles,
            "header_pstyles": ParagraphStyle("cols", fontSize=9),
            "cell_styles": cell_styles,
            "header_styles": header_styles,
            "col_widths": list(Widths.sv_snp.values()),
        }
        # Report is in the following order
        # TODO: create universal footer/header spec for tables
        # TODO: build the pdf for Therapeutic information
        # TODO: build the pdf for Copy number
        # TODO: build the pdf for signature analysis
        # TODO: create universal styles for variant tables
        # front_page: ReportElement =
        table_decorator = None
        self.build_table(
            table_spec,
            vtable_styles,
            self.data["relevant"]["small"],
            "Relevant Small Variants",
            "Relevant Small Variants (Continued)",
            table_decorator,
            "rs.pdf",
        )
        self.build_table(
            table_spec,
            vtable_styles,
            self.data["nonrelevant"]["small"],
            "Non-relevant Small Variants",
            "Non-relevant Small Variants (Continued)",
            table_decorator,
            "nrs.pdf",
        )
        self.build_table(
            table_spec,
            vtable_styles,
            self.data["nonrelevant"]["sv"],
            "Non-relevant Structural Variants",
            "Non-relevant Structural Variants (Continued)",
            table_decorator,
            "rv.pdf",
        )
        self.build_table(
            table_spec,
            vtable_styles,
            self.data["nonrelevant"]["sv"],
            "Non-relevant Structural Variants",
            "Non-relevant Structural Variants (Continued)",
            table_decorator,
            "nrv.pdf",
        )


class ReportElement:
    """Generic class to use to add report elements
    These are intended to be discrete units of the report e.g. a specific table or
    front page

    Used to give more control over how flowables and canvas interact with one another
    e.g. making a table that has repeated headers that change depending on the page number

    Important fields in `spec`
    - [left|right|bottom|top]_margin
    - header|footer: text for the header and footer
    - [header|footer]_[first|later]: this overrides header and footer
    - [header|footer]_pos: positions of header and footer
    - [header|footer]_style: ParagraphStyle used for headers and footers
    """

    def __init__(
        self,
        filename: str,
        spec: dict,
        pagesize=A4,
        w_page_break=True,
    ) -> None:
        self.pdf = SimpleDocTemplate(
            filename,
            pagesize=pagesize,
            rightMargin=spec.get("right_margin", inch),
            leftMargin=spec.get("left_margin", inch),
            topMargin=spec.get("top_margin", inch),
            bmargin=spec.get("bottom_margin", inch),
        )
        self.spec = spec
        self.width = pagesize[0]
        self.decorator: Callable[[Canvas, BaseDocTemplate], None] = None
        self.height = pagesize[1]
        self.w_page_break = w_page_break
        self.elements = []

    def add_decorator(self, fn: Callable[[Canvas, BaseDocTemplate], None]) -> None:
        self.decorator = fn

    def draw_pages(self, h, f) -> Callable:
        """Generic function for drawing pages with header, footer if they are provided,
        as well as instance-specific decorations

        :param h: header text (can also be provided in spec dict)
        :param f: footer text (can also be provided in spec dict)

        :returns: an anonymous function to use by ReportLab to build the pdf
        """

        def fn(c, d):
            if self.decorator:
                self.decorator(c, d)
            if htext := self.spec.get("header", h):
                draw_paragraph(
                    htext,
                    self.spec.get(
                        "header_pos", (self.width - 5 * inch, self.height - inch)
                    ),
                    self.spec.get("header_style", STYLES["h3"]),
                    c,
                    d,
                )
            # Draw the footer
            if ftext := self.spec.get("footer", f):
                draw_paragraph(
                    ftext,
                    self.spec.get("footer_pos", (self.width - 5 * inch, 0 + inch)),
                    self.spec.get("footer_style", STYLES["Normal"]),
                    c,
                    d,
                )

        return fn

    def add_table(
        self,
        data: pl.DataFrame,
        cell_pstyles: ParagraphStyle | dict[int, ParagraphStyle],
        header_pstyles: ParagraphStyle | dict[int, ParagraphStyle],
        cell_styles: list,
        header_styles: list,
        col_widths,
        repeat_header=True,
    ) -> None:
        """Add a LongTable to elements, formatting accordingly

        :param data_styles: Dictionary or ParagraphStyle specifying styles to apply to
        each column. See signature of fn  `add_pstyles`
        :param header_styles: Dictionary specifying styles to apply to only to the
        field (column) headers

        """
        cells = add_pstyles(data, cell_pstyles)
        headers = add_pstyles(data.columns, header_pstyles)
        cells.insert(0, headers)
        table: LongTable = LongTable(
            cells,
            colWidths=col_widths,
            repeatRows=1 if repeat_header else 0,
        )
        table.setStyle(cell_styles + header_styles)
        self.elements.append(table)

    def build(self) -> None:
        if self.w_page_break:
            self.elements.append(PageBreak())
        self.pdf.build(
            self.elements,
            onFirstPage=self.draw_pages(
                self.spec.get("header_first", ""),
                self.spec.get("footer_first", ""),
            ),
            onLaterPages=self.draw_pages(
                self.spec.get("header_later", ""),
                self.spec.get("footer_later", ""),
            ),
        )
