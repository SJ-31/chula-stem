import os
from typing import Callable

import polars as pl
from pymupdf import Document
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

from chula_stem.report.spec import BOLD_FONT

STYLES = getSampleStyleSheet()


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


# * Utility fns
#


def no_data_decorator(c: Canvas, d: Document) -> None:
    style = ParagraphStyle(
        "nd", fontName=BOLD_FONT, fontSize=40, borderColor=colors.black
    )
    # TODO: try to contain the text in a rectangle
    # c.saveState()
    # c.rect(A4[0] / 2, A4[1] / 2, 5 * cm, 4 * cm)
    # c.restoreState()
    draw_paragraph("NO DATA", (A4[0] / 2 - 3 * cm, A4[1] / 2 + 2 * inch), style, c, d)


def draw_paragraph(
    text: str,
    pos: tuple,
    style: ParagraphStyle,
    canvas: Canvas,
    doc: BaseDocTemplate,
) -> None:
    canvas.saveState()
    para = Paragraph(text, style)
    para.wrap(doc.width, doc.topMargin)
    para.drawOn(canvas, pos[0], pos[1])
    canvas.restoreState()


def add_pstyles(
    data: pl.DataFrame | list,
    style: ParagraphStyle | dict[int, ParagraphStyle],
    nested=False,
) -> list:
    """Helper for adding paragraph styles to lists or data in dfs

    :param style: A single style which is then applied to all data.
    Alternatively, a map of column_index > style specifying styles to
    apply to specific columns. A key for 'None' is the default and applied to columns
    not explicitly given

    :returns:
    """
    if isinstance(style, dict):
        if not style.get(None):
            raise ValueError("A default key `None` must be provided!")

        style_fn = lambda x, index: Paragraph(str(x), style.get(index, style[None]))
    else:
        style_fn = lambda x, _: Paragraph(str(x), style)

    if isinstance(data, pl.DataFrame):
        return [[style_fn(s, i) for i, s in enumerate(row)] for row in data.iter_rows()]
    elif nested:
        return [[style_fn(s, i) for i, s in enumerate(row)] for row in data]
    else:
        return [style_fn(row, i) for i, row in enumerate(data)]


class ResultsReport:
    def __init__(
        self, filename: str, font: str, bold_font: str, tmpdir: str = "temp"
    ) -> None:
        self.filename: str = filename
        self.tmpdir: str = tmpdir
        self.font: str = font
        self.bold_font: str = bold_font
        try:
            os.makedirs(self.tmpdir)
        except FileExistsError:
            print("WARNING: directory exists")
        os.chdir(self.tmpdir)

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


# * Report element
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
        self.decorators: list[Callable[[Canvas, BaseDocTemplate], None]] = []
        self.height = pagesize[1]
        self.w_page_break = w_page_break
        self.elements = []

    def add_decorator(self, fn: Callable[[Canvas, BaseDocTemplate], None]) -> None:
        self.decorators.append(fn)

    def draw_pages(self, h, f) -> Callable:
        """Generic function for drawing pages with header, footer if they are provided,
        as well as instance-specific decorations

        :param h: header text (can also be provided in spec dict)
        :param f: footer text (can also be provided in spec dict)

        :returns: an anonymous function to use by ReportLab to build the pdf
        """

        def fn(c, d):
            for decorator in self.decorators:
                if decorator:
                    decorator(c, d)
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
        if header_pstyles and header_styles:
            headers = add_pstyles(data.columns, header_pstyles)
            cells.insert(0, headers)
        table: LongTable = LongTable(
            cells,
            colWidths=col_widths,
            repeatRows=1 if repeat_header else 0,
        )
        style_list = []
        if cell_styles:
            style_list.extend(cell_styles)
        if header_styles:
            style_list.extend(header_styles)
        if style_list:
            table.setStyle(style_list)
        if data.is_empty():
            self.add_decorator(no_data_decorator)
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
