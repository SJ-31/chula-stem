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
    VTABLE_COL_WIDTHS,
    CTABLE_WIDTHS,
    VTABLE_COL_WIDTHS,
    CTABLE_RENAME,
    VTABLE_RENAME,
    URLS,
)

STYLES = getSampleStyleSheet()

from chula_stem.callset_qc import IMPACT_MAP
from chula_stem.databases import add_therapy_info
from chula_stem.utils import add_loc, read_facets_rds


def add_link(text: str, link: str, **kwargs) -> str:
    """
    Wrap text in RML link tag. Can use any tag defined in
        https://docs.reportlab.com/tagref/tag_link/
    """
    begin = f'<link href="{link}"'
    if kwargs:
        for k, v in kwargs.items():
            begin = f'{begin} {k}="{v}"'
    begin = f"{begin}>"
    return f"{begin}{text}</link>"


def get_copy_number(
    df: pl.DataFrame, facets_path: str = "", cnvkit_path: str = "", cn_col_name="cn"
) -> pl.DataFrame:
    """
    Add a copy number column `cn_col_name` to df
    by merging with the raw CNV data. df must have the
    column "Locus", containing the location of the cnv in the form of
        chromosome:start-end
    :param format: one of CNVKit|Facets
    """
    wanted_cols = df.columns + [cn_col_name]
    if facets_path:
        facets: pl.DataFrame = read_facets_rds(facets_path).rename(
            {"tcn.em": cn_col_name}
        )
        facets = add_loc(facets, "chrom")
        df = df.join(facets, how="left", on="Locus")
    if cnvkit_path:
        cnvkit: pl.DataFrame = pl.read_csv(
            cnvkit_path, separator="\t", null_values="NA"
        )
        cnvkit = add_loc(cnvkit)
        df = df.join(cnvkit, how="left", on="Locus")
    if facets_path and cnvkit_path:
        df = df.with_columns(
            pl.coalesce([cn_col_name, f"{cn_col_name}_right"]).alias(cn_col_name)
        )
    return df.select(wanted_cols)


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
        call_all: bool = True,
    ) -> None:
        self.civic_cache = civic_cache
        self.pandrugs2_cache = pandrugs2_cache
        self.data: dict = {"relevant": {}, "nonrelevant": {}, "all": {}}
        self.tmpdir = tmpdir

        if call_all:
            try:
                os.makedirs(self.tmpdir)
            except FileExistsError:
                print("WARNING: directory exists")
            os.chdir(self.tmpdir)
            self._format_vep(vep_small, "small")
            self._format_vep(vep_sv, "sv")
            self._format_classifycnv(classify_cnv, facets, cnvkit)

    def _format_classifycnv(
        self,
        classifycnv_path: str,
        facets_path: str = "",
        cnvkit_path: str = "",
    ) -> None:
        cnv = pl.read_csv(
            classifycnv_path, separator="\t", null_values="NA", infer_schema_length=None
        ).with_columns(
            pl.struct(s="source", acc="accession")
            .map_elements(
                lambda x: (
                    self._format_dbvar_link(x["acc"]) if x["acc"] != "NA" else x["s"]
                ),
                return_dtype=pl.String,
            )
            .alias("source")
        )
        wanted_cols = ["Locus"] + list(CTABLE_RENAME.values())
        cnv = (
            cnv.with_columns(
                pl.concat_str(["Start", "End"], separator="-").alias("range"),
            )
            .with_columns(
                pl.concat_str(["Chromosome", "range"], separator=":").alias("Locus"),
            )
            .rename(CTABLE_RENAME)
        ).select(wanted_cols)
        cn_col = "Estimated Copy Number"
        wanted_cols.insert(1, cn_col)
        with_cn = get_copy_number(
            cnv, facets_path, cnvkit_path, cn_col_name="cn"
        ).rename({"cn": cn_col}).selekt(wanted_cols)
        relevant = with_cn.filter(
            pl.col("Known/predicted dosage-sensitive genes").is_not_null()
        )
        nonrelevant = with_cn.filter(
            pl.col("Known/predicted dosage-sensitive genes").is_null()
        )
        self.data["relevant"]["cnv"] = relevant
        self.data["nonrelevant"]["cnv"] = nonrelevant
        self.data["all"]["cnv"] = with_cn

    @staticmethod
    def _format_dbvar_link(x) -> str:
        link = add_link(x, f"{URLS['dbvar']}/{x}", underline="yes")
        return f"dbVar:{link}"

    def _format_vep(self, vep_path: str, variant_class: str) -> None:
        """Format and filter vep output into a dataframe with values ready to write
        into a reportlab table
        """
        r_out = f"{self.tmpdir}/{variant_class}_relevant.parquet"
        nr_out = f"{self.tmpdir}/{variant_class}_non-relevant.parquet"
        all_out = f"{self.tmpdir}/{variant_class}_all.parquet"
        if os.path.exists(r_out) and os.path.exists(nr_out) and os.path.exists(all_out):
            # self._read_vep_saved(r_out, nr_out, all_out, variant_class)
            for c, p in zip(
                ["relevant", "nonrelevant", "all"], [r_out, nr_out, all_out]
            ):
                self.data[c][variant_class] = pl.read_parquet(p)
            return
        var_col: str = "Database Name"  # New column for 'Existing_variation'
        with_therapeutics: pl.DataFrame = add_therapy_info(
            vep_path, self.civic_cache, self.pandrugs2_cache
        )
        wanted_cols: list = list(VTABLE_RENAME.values())
        with_therapeutics = (
            with_therapeutics.drop("Gene")
            .rename(VTABLE_RENAME)
            .with_columns(
                pl.col("HGVS").str.extract(r".*:(.*)$", 1),
                pl.col("Variant Type").str.replace("_variant$", ""),
                (pl.col("Variant Allele Frequency").cast(pl.Float64) * 100)
                .round(3)
                .map_elements(lambda x: f"{x}%", return_dtype=pl.String),
            )
            .with_columns(
                cs.by_dtype(pl.String)
                .fill_null("-")
                .str.replace_all("&", ", ", literal=True)
                .str.replace_all("_", " ", literal=True),
            )
        )
        confident = (
            with_therapeutics.filter(
                (pl.col(var_col).is_not_null()) & (pl.col("SOMATIC") == "1")
            )
            .with_columns(impact_score=pl.col("IMPACT").replace_strict(IMPACT_MAP))
            .sort(pl.col("impact_score"), descending=True)
        )
        others = with_therapeutics.filter(~pl.col("VAR_ID").is_in(confident["VAR_ID"]))
        self.data["relevant"][variant_class] = confident.select(wanted_cols)
        self.data["nonrelevant"][variant_class] = others.select(wanted_cols)
        self.data["all"][variant_class] = with_therapeutics
        for category, path in zip(
            ["relevant", "nonrelevant", "all"], [r_out, nr_out, all_out]
        ):
            self.data[category][variant_class].write_parquet(path)

            # self.data[category][variant_class].with_columns(
            #     cs.by_dtype(pl.List(pl.String)).list.join(";")
            # ).write_csv(path, separator="\t", null_value="NA")

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
            "col_widths": list(VTABLE_COL_WIDTHS.values()),
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
