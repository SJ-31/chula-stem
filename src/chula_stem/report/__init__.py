import os
from typing import Callable

import polars as pl
import polars.selectors as cs
from pymupdf import Document, Page
import pymupdf
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
from chula_stem.databases import get_therapy_df
from chula_stem.report.spec import (
    cnv_style,
    reference_list_style,
    signature_style,
    sv_snp_style,
    repeat_style,
    therapy_style,
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
    data: pl.DataFrame | list, style: ParagraphStyle | dict[int, ParagraphStyle]
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


class VariantCallingReport(ResultsReport):
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
        msisensor_pro: str = "",
        sigprofiler: str = "",
        cosmic_reference: str = "",
        tmpdir: str = "temp",
        font="Helvetica",
        bold_font="Helvetica-Bold",
    ) -> None:
        super().__init__(
            filename=filename,
            tmpdir=tmpdir,
            font=font,
            bold_font=bold_font,
        )
        self.civic_cache: str = civic_cache
        self.header_y: int = cm
        self.pandrugs2_cache: str = pandrugs2_cache
        self.data: dict = {"relevant": {}, "nonrelevant": {}, "all": {}}
        self.files: dict = {"signatures": []}
        self.table_styles: dict = {
            "small": sv_snp_style(),
            "cnv": cnv_style(),
            "repeat": repeat_style(),
        }
        self.table_styles["sv"] = self.table_styles["small"]
        # TODO: add in the style for therapies

        gene_spec = [
            {"file": vep_small, "column": "SYMBOL"},
            {"file": vep_sv, "column": "SYMBOL"},
            {
                "file": classify_cnv,
                "is_list": True,
                "is_list_separator": ",",
                "column": "All protein coding genes",
            },
            {"file": msisensor_pro, "column": "gene_name"},
        ]
        therapy_data_cache: str = "therapy_data.parquet"
        if os.path.exists(therapy_data_cache):
            therapy_df = pl.read_parquet(therapy_data_cache)
        else:
            gene_list = fr.get_genes(gene_spec)
            therapy_df: pl.DataFrame = get_therapy_df(
                gene_list, civic_cache=civic_cache, pandrugs2_cache=pandrugs2_cache
            )
            therapy_df.write_parquet(therapy_data_cache)

        calls = [
            lambda: fr.vep_fmt(vep_small, tmpdir, "small", therapy_df),
            lambda: fr.vep_fmt(vep_sv, tmpdir, "sv", therapy_df),
            lambda: fr.classify_cnv_fmt(classify_cnv, facets, cnvkit),
            lambda: fr.msisensor_pro_fmt(msisensor_pro),
        ]
        self.variant_names: dict[str, str] = dict(
            zip(
                ["small", "sv", "cnv", "repeat"],
                [
                    "Small Variants",
                    "Structural Variants",
                    "Copy Number Variants",
                    "Tandem Repeats",
                ],
            )
        )

        variant_spec: list[dict] = []
        for type, fn_call in zip(self.variant_names, calls):
            all, relevant, nonrelevant = fn_call()
            self.data["relevant"][type] = relevant
            self.data["nonrelevant"][type] = nonrelevant
            self.data["all"][type] = all
            spec: dict = {
                "type": self.variant_names[type],
                "df": relevant,
                "is_list": False,
                "gene_col": "Gene",
            }
            if type == "cnv":
                spec["is_list"] = True
                spec["gene_col"] = "All Genes"
                spec["separator"] = ","
            elif type == "repeat":
                spec["gene_col"] = "Affected Gene"
            variant_spec.append(spec)

        self.data["therapies"], self.data["study_references"] = fr.therapy_fmt(
            therapy_df, variant_spec
        )
        self.data["sigprofiler"] = fr.sigprofiler_fmt(sigprofiler, cosmic_reference)

    def build(self) -> None:
        """Create pdf files for all report elements individually, concatenate them
        and add header (time + page number)
        """
        table_spec = {"header_pos": (A4[0] - 5 * inch, A4[1] - inch)}
        # Report is in the following order
        # TODO: create universal footer/header spec for tables
        # TODO: create universal styles for variant tables
        # front_page: ReportElement =
        table_decorator = None
        for table, name in self.variant_names.items():
            style = self.table_styles[table]
            self.build_table(
                table_spec,
                style,
                self.data["relevant"][table],
                f"Relevant {name}",
                f"Relevant {name} (Continued)",
                table_decorator,
                f"rel_{table}.pdf",
            )
            self.build_table(
                table_spec,
                style,
                self.data["nonrelevant"][table],
                f"Non-relevant {name}",
                f"Non-relevant {name} (Continued)",
                table_decorator,
                f"non-rel_{table}.pdf",
            )
        self.build_table(
            table_spec,
            therapy_style(),
            self.data["therapies"],
            "Relevant therapies",
            "Relevant therapies (Continued)",
            table_decorator,
            "therapies.pdf",
        )
        self.build_table(
            table_spec,
            reference_list_style(),
            self.data["study_references"],
            "References",
            "",
            table_decorator,
            "reference_list.pdf",
        )
        for n, sample in enumerate(self.data["sigprofiler"]):
            if len(self.data["sigprofiler"]) == 1:
                first = "Mutational signatures"
                last = "Mutational signatures (continued)"
            else:
                first = f"Mutational signatures, Sample {n}"
                last = f"Mutational signatures, Sample {n} (continued)"
            name = f"signatures_{n}.pdf"
            self.build_table(
                table_spec,
                signature_style(),
                sample,
                first,
                last,
                table_decorator,
                name,
            )
            self.files["signatures"].append(name)

    def _add_page_numbers(self, page: Page, total_pages: int, offset: int = 1) -> None:
        n: int = page.number
        text = f"Page {n + offset} of {total_pages}"
        page.insert_text(
            (A4[0] - inch, A4[1] - inch),
            text,
            fontsize=20,
            fontname="Courier-Bold",
        )

    # TODO: add in the front, back pages
    def merge(self) -> None:
        toc: list = []
        doc: Document = pymupdf.open()

        toc.append([1, "Relevant therapies", len(doc) + 1])
        doc.insert_file("therapies.pdf")

        for t, n in self.variant_names.items():
            toc.append([1, n, len(doc) + 1])
            toc.append([2, "Relevant", len(doc) + 1])
            doc.insert_file(f"rel_{t}.pdf")
            toc.append([2, "Non-relevant", len(doc) + 1])
            doc.insert_file(f"non-rel_{t}.pdf")

        toc.append([1, "Mutational signatures", len(doc) + 1])
        if len(self.files["signatures"]) == 1:
            doc.insert_file(self.files["signatures"][0])
        else:
            for index, s in enumerate(self.files["signatures"]):
                toc.append([2, f"Sample {index}", len(doc) + 1])
                doc.insert_file(s)

        toc.append([1, "References", len(doc) + 1])
        doc.insert_file("reference_list.pdf")

        total_pages: int = len(doc)
        for p in doc.pages():
            self._add_page_numbers(p, total_pages)
        doc.set_toc(toc)
        doc.save(self.filename)


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
