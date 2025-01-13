#!/usr/bin/env ipython
import datetime
import os
from pathlib import Path

import polars as pl
import pymupdf
from pymupdf import Document, Page
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import cm, inch
from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus import (
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
)

import chula_stem.report.format as fr
from chula_stem.databases import get_therapy_df
from chula_stem.plotting import plot_cnvkit
from chula_stem.report import (
    PATIENT_DATA,
    SAMPLE_DATA,
    ResultsReport,
    add_pstyles,
    draw_paragraph,
)
from chula_stem.report.spec import (
    AVAILABLE_WIDTH,
    BOLD_FONT,
    FONT,
    STYLE,
    Rename,
    cnv_style,
    reference_list_style,
    repeat_style,
    signature_style,
    snp_style,
    sv_style,
    therapy_style,
)
from chula_stem.report.utils import style_cells

HGVS_LINK: str = "https://hgvs-nomenclature.org/stable"

# * Add acronyms
ACRONYMS = [
    (
        "<b>• VRS:</b>",
        "Variant Read Support (raw number of reads mapping to alternate allele)",
    ),
    ("<b>• VAF:</b>", "Variant Allele Frequency"),
    (
        "<b>• HGVS:</b>",
        fr.add_link("Human Genome Variation Society nomenclature", HGVS_LINK),
    ),
]
ACRONYM_STYLE = ParagraphStyle("acronyms", fontSize=8, fontName=FONT)
ACRONYM_TABLE: Table = Table(
    add_pstyles(ACRONYMS, ACRONYM_STYLE, True),
    style=style_cells((0, 0), valign="TOP"),
    colWidths=(AVAILABLE_WIDTH * 0.09, AVAILABLE_WIDTH * 0.91),
    rowHeights=11,
)


def override_style(style: ParagraphStyle, **kwargs) -> ParagraphStyle:
    """Return a new ParagraphStyle, inheriting attributes from `style`
    and overriding any attributes specified in `override`
    """
    inherit: dict = {k: getattr(style, k) for k in dir(style) if not k.startswith("__")}
    inherit.update(kwargs)
    return ParagraphStyle(**inherit)


# * Add criteria info
CRITERIA_STYLE = ParagraphStyle("criteria", fontSize=8, fontName=FONT)
CRITERIA = [
    (
        "",
        "<b>In 'Relevant' table</b>",
        "<b>All entries</b>",
    ),
    (
        "<b>Small and Structural Variants</b>",
        'Must have a ClinVar annotation to "pathogenic" or "likely pathogenic"',
        'Removed entries having only ClinVar annotations of "benign" or "likely benign"',
    ),
    (
        "<b>Copy Number Variants</b>",
        "At least one gene in the interval must have ClinGen dosage-sensitive information",
        "None",
    ),
    (
        "<b>Tandem repeats</b>",
        "Affected gene must have ClinGen disease-association information",
        "Repeats must overlap known genes",
    ),
]
CRITERIA_CELLS = (
    style_cells((0, 0), valign="TOP")
    + style_cells((0, 1), end=(2, 1), background="#c3e6e9")
    + style_cells((0, 3), end=(2, 3), background="#c3e6e9")
)
CRITERIA_TABLE: Table = Table(
    add_pstyles(CRITERIA, CRITERIA_STYLE, True),
    style=CRITERIA_CELLS,
    colWidths=(AVAILABLE_WIDTH * 0.3, AVAILABLE_WIDTH * 0.4, AVAILABLE_WIDTH * 0.3),
)


class VariantCallingReport(ResultsReport):
    def __init__(
        self,
        filename: str,
        metadata: dict,
        other_text: dict = {},
        pandrugs2_cache: str = "pandrugs2_cache.json",
        civic_cache: str = "civic_cache.json",
        therapy_data_cache: str = "therapy_data.parquet",
        vep_small: str = "",
        vep_sv: str = "",
        classify_cnv: str = "",
        facets: str = "",
        cnvkit_cns: str = "",
        cnvkit_cnr: str = "",
        msisensor_pro: str = "",
        sigprofiler: str = "",
        cosmic_reference: str = "",
        tmpdir: str = "temp",
        font=FONT,
        bold_font=BOLD_FONT,
        plot=True,
        filter_confident_therapies=True,
        show_therapies=True,
    ) -> None:
        super().__init__(
            filename=filename,
            metadata=metadata,
            tmpdir=tmpdir,
            font=font,
            bold_font=bold_font,
            other_text=other_text,
        )
        self.spec: dict = {
            "right_margin": 2 * cm,
            "left_margin": 2 * cm,
            "top_margin": 2 * cm,
            "bottom_margin": 2 * cm,
            "header_pos": (2 * cm, A4[1] - 1.5 * cm),
            "footer_pos": (A4[0] - 2 * cm, A4[1] - cm),
            "header_style": STYLE["table_title_style"],
        }
        self.show_therapies = show_therapies
        self.civic_cache: str = civic_cache
        self.pandrugs2_cache: str = pandrugs2_cache
        self.data: dict = {"relevant": {}, "nonrelevant": {}, "all": {}}
        self.files: dict = {"signatures": []}
        self.table_styles: dict = {
            "small": snp_style(),
            "cnv": cnv_style(),
            "repeat": repeat_style(),
            "sv": sv_style(),
        }
        if plot:
            self.plot_cnvkit(cnvkit_cnr, cnvkit_cns)

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
        if os.path.exists(therapy_data_cache) and self.show_therapies:
            therapy_df = pl.read_parquet(therapy_data_cache)
        elif self.show_therapies:
            gene_list = fr.get_genes(gene_spec)
            therapy_df: pl.DataFrame = get_therapy_df(
                gene_list,
                civic_cache=civic_cache,
                pandrugs2_cache=pandrugs2_cache,
                filter_confident=filter_confident_therapies,
            )
            therapy_df.write_parquet(therapy_data_cache)
        else:
            therapy_df = pl.DataFrame()

        calls = [
            lambda: fr.vep_fmt(vep_small, "small"),
            lambda: fr.vep_fmt(vep_sv, "sv"),
            lambda: fr.classify_cnv_fmt(classify_cnv, facets, cnvkit_cns),
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
        self.stats: dict = {"counts": {}}
        for type, fn_call in zip(self.variant_names, calls):
            all, relevant, nonrelevant = fn_call()
            self.data["relevant"][type] = relevant
            formatted_name = self.variant_names.get(type)
            self.stats["counts"][f"Relevant {formatted_name}"] = relevant.shape[0]
            self.stats["counts"][f"Non-relevant {formatted_name}"] = nonrelevant.shape[
                0
            ]
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
                spec["gene_col"] = Rename.cnv["All protein coding genes"]
                spec["separator"] = ","
            elif type == "repeat":
                spec["gene_col"] = Rename.repeat["gene_name"]
            variant_spec.append(spec)

        if self.show_therapies:
            self.data["therapies"], self.data["study_references"] = fr.therapy_fmt(
                therapy_df, variant_spec
            )
            self.stats["counts"]["Available Therapies"] = self.data["therapies"].shape[
                0
            ]
        self.data["sigprofiler"] = fr.sigprofiler_fmt(sigprofiler, cosmic_reference)
        self.stats["counts"]["Mutational Signatures"] = sum(
            [d.shape[0] for d in self.data["sigprofiler"]]
        )

    def build_front_page(self) -> None:
        """Front page
        Shows the following:
        - Patient information: name, date of birth, gender, hospital, diagnosis
        - Sample information: id, collection date, sample type, cancer type
        - Summary
        """

        def decorator(c: Canvas, d: SimpleDocTemplate) -> None:
            c.saveState()
            # Add decorations
            c.setFillColor("#86cfd5")
            rect_y = A4[1] - cm * 4.5
            c.rect(0, rect_y, A4[0], cm * 4.5, stroke=0, fill=1)
            c.setLineWidth(0.9)
            c.setStrokeColor("#86cfd5")
            y1 = cm * 5
            y2 = A4[1] - cm * 10
            x, xend = cm * 2, A4[0] - cm * 2
            c.line(x, y1, xend, y1)
            c.line(x, y2, xend, y2)
            c.restoreState()

            note = self.get_other_text("front_page_note", "")
            disclaimer = self.get_other_text("disclaimer", "")
            title = self.get_other_text("title", "")
            details = self.get_other_text("front_page_details", "")

            # Notes and disclaimer
            note_style = ParagraphStyle("note", fontSize=7, fontName=FONT)
            note_pos = (x, cm * 2 + 10)
            draw_paragraph(f"<b>Note:</b><br/>{note}", note_pos, note_style, c, d)
            disclaimer_pos = (x, cm)
            draw_paragraph(
                f"<b>Disclaimer:</b> {disclaimer}", disclaimer_pos, note_style, c, d
            )

            # Title
            title_y = A4[1] - 3 * cm
            draw_paragraph(
                title,
                (x, title_y),
                STYLE["title_style"],
                c,
                d,
                width=AVAILABLE_WIDTH * 0.4,
            )
            draw_paragraph(
                details,
                (2.8 * cm + AVAILABLE_WIDTH * 0.4, title_y),
                STYLE["detail_style"],
                c,
                d,
                width=AVAILABLE_WIDTH * 0.6,
            )

        # Patient, sample info tables
        fontsize = 8
        header_fontsize = 13
        info_table_widths = [100, 180]
        cell_row_height = 15
        styles = {
            # 0: ParagraphStyle("keys", fontName=BOLD_FONT, fontSize=fontsize),
            None: ParagraphStyle("values", fontName=FONT, fontSize=fontsize),
        }
        patient_data, sample_data = [], []
        for spec, add_to in zip(
            [PATIENT_DATA, SAMPLE_DATA],
            [patient_data, sample_data],
        ):
            for k, v in spec.items():
                if k not in self.ignored_data_fields:
                    row = [f"{k}:", self.metadata.get(v, "-")]
                    add_to.append(row)

        patient_data = add_pstyles(patient_data, styles, nested=True)
        sample_data = add_pstyles(sample_data, styles, nested=True)
        cell_styles = style_cells((0, 0), ncols=2, valign="TOP")
        header_pstyle = ParagraphStyle(
            "header",
            fontName=BOLD_FONT,
            fontSize=header_fontsize,
            textColor="#86cfd5",
            firstLineIndent=7,
        )
        p_header, s_header = (
            add_pstyles(["Patient Information"], header_pstyle),
            add_pstyles(["Sample Information"], header_pstyle),
        )

        patient_table = Table(
            patient_data, colWidths=info_table_widths, rowHeights=cell_row_height
        )
        patient_table.setStyle(cell_styles)
        sample_table = Table(
            sample_data, colWidths=info_table_widths, rowHeights=cell_row_height
        )
        sample_table.setStyle(cell_styles)
        spacing = Spacer(AVAILABLE_WIDTH, cm)
        stacked = Table(
            [p_header, [patient_table], [spacing], s_header, [sample_table]],
        )
        full = Table(
            [
                [Spacer(AVAILABLE_WIDTH, 3 * cm)],
                ["", stacked],
            ],
            colWidths=[AVAILABLE_WIDTH * 0.4, AVAILABLE_WIDTH * 0.6],
        )

        doc = SimpleDocTemplate(
            "front_page.pdf",
            pagesize=A4,
            rightMargin=self.spec.get("right_margin", inch),
            leftMargin=self.spec.get("left_margin", inch),
            topMargin=self.spec.get("top_margin", inch),
            bmargin=self.spec.get("bottom_margin", inch),
        )
        doc.build(
            [full],
            onFirstPage=decorator,
        )

    def build_toc(self, toc: list) -> None:
        """Build a toc page

        :param toc: toc list of the form used by pymupdf
        """
        cell_sets = [[], []]
        cell_styles = []
        cur_set = cell_sets[0]
        i, offset = 0, 0
        for level, text, page in toc:
            t = Paragraph(text, STYLE["toc_text"](level))
            n = Paragraph(str(page), STYLE["toc_num"])

            tstyle = style_cells(
                (offset, i), end=(offset + 2, i), background=STYLE["toc_cell_bg"]
            )
            if i % 2 == 0:
                cell_styles.extend(tstyle)
            i += 1
            if text == "Tandem Repeats":
                i = 0
                offset = 2
                cur_set = cell_sets[1]
            cur_set.append([t, n])
        cell_styles.extend(style_cells((2, 0), end=(2, -1), background=colors.white))
        n = max(len(cell_sets[0]), len(cell_sets[1]))
        cell_sets[0] += [""] * (n - len(cell_sets[0]))
        cell_sets[1] += [""] * (n - len(cell_sets[1]))

        all_cells = []
        for left, right in zip(*cell_sets):
            if isinstance(left, list) and isinstance(right, list):
                all_cells.append([*left, " ", *right])
            elif isinstance(left, list):
                all_cells.append([*left, "", "", ""])
            else:
                all_cells.append(["", "", "", *right])

        widths = [((AVAILABLE_WIDTH / 4) - 3)] * 4
        widths.insert(2, 12)
        toc_table = Table(all_cells, colWidths=widths)
        toc_table.setStyle(style_cells((0, 0), valign="TOP") + cell_styles)

        doc = SimpleDocTemplate(
            "toc.pdf",
            pagesize=A4,
            rightMargin=self.spec.get("right_margin", inch),
            leftMargin=self.spec.get("left_margin", inch),
            topMargin=self.spec.get("top_margin", inch),
            bmargin=self.spec.get("bottom_margin", inch),
        )
        stats = [[v, "", k] for k, v in self.stats["counts"].items()]
        stats_table = Table(stats)
        # full = Table([[toc_table, "", stats_table]], colWidths=[])

        def decorator(c: Canvas, d: SimpleDocTemplate) -> None:
            STYLE["table_decorator"](c, d)
            draw_paragraph(
                "Table of Contents",
                self.spec.get("header_pos"),
                STYLE["toc_title_style"],
                c,
                d,
            )
            c.saveState()
            c.setLineWidth(0.9)
            c.setStrokeColor("#86cfd5")
            y1 = cm * 21
            c.line(cm * 2, y1, cm * 2 + AVAILABLE_WIDTH, y1)

        doc.build(
            [
                toc_table,
                Spacer(width=0, height=1.6 * cm),
                Paragraph(
                    "Acronyms",
                    override_style(ACRONYM_STYLE, fontSize=12, textColor="#549e9e"),
                ),
                Spacer(width=0, height=5),
                ACRONYM_TABLE,
                Spacer(width=0, height=1.5 * cm),
                Paragraph(
                    "Criteria for variant inclusion",
                    override_style(CRITERIA_STYLE, fontSize=12, textColor="#549e9e"),
                ),
                Spacer(width=0, height=5),
                CRITERIA_TABLE,
            ],
            onFirstPage=decorator,
        )

    @staticmethod
    def plot_cnvkit(cnr: str, cns: str) -> None:
        sizes = {"ncol": 3, "height": 30, "dpi": 200, "width": 20}
        file = "cnvkit.png"
        plot_cnvkit(cnr, cns, None, sizes, file)
        doc: Document = pymupdf.open()
        page = doc.new_page()
        with open(file, "rb") as pic:
            img = pic.read()
            page.insert_image(page.rect, stream=img)
        doc.save("cnvkit.pdf")

    def build(self) -> None:
        """Create pdf files for all report elements individually, concatenate them
        and add header (time + page number)
        """
        table_decorator = STYLE["table_decorator"]
        self.build_front_page()

        for table, name in self.variant_names.items():
            style = self.table_styles[table]
            self.build_table(
                self.spec,
                style,
                self.data["relevant"][table],
                f"Relevant {name}",
                f"Relevant {name} (Continued)",
                table_decorator,
                f"rel_{table}.pdf",
            )
            self.build_table(
                self.spec,
                style,
                self.data["nonrelevant"][table],
                f"Non-relevant {name}",
                f"Non-relevant {name} (Continued)",
                table_decorator,
                f"non-rel_{table}.pdf",
            )
        if self.show_therapies:
            self.build_table(
                self.spec,
                therapy_style(),
                self.data["therapies"],
                "Relevant therapies",
                "Relevant therapies (Continued)",
                table_decorator,
                "therapies.pdf",
            )
            if not self.data["study_references"].is_empty():
                self.build_table(
                    self.spec,
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
                self.spec,
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
        number = f"| {n + offset} of {total_pages}"
        date = datetime.datetime.now().strftime(
            "%Y-%m-%d",
        )
        x, y = self.spec["footer_pos"]
        size = 8
        page.insert_text(
            (x, y),
            number,
            fontsize=size,
            fontname=FONT,
        )
        page.insert_text(
            (x - 1.6 * cm, y),
            date,
            fontsize=size,
            fontname=FONT,
            color=(0.34902, 0.63529, 0.63529),
        )

    def merge(self) -> None:
        toc: list = []
        doc: Document = pymupdf.open()

        offset = 2  # Account for front page and TOC page
        doc.insert_file("front_page.pdf")
        if self.show_therapies:
            toc.append([1, "Relevant therapies", len(doc) + offset])
            doc.insert_file("therapies.pdf")

        for t, n in self.variant_names.items():
            toc.append([1, n, len(doc) + offset])
            toc.append([2, "Relevant", len(doc) + offset])
            doc.insert_file(f"rel_{t}.pdf")
            toc.append([2, "Non-relevant", len(doc) + offset])
            doc.insert_file(f"non-rel_{t}.pdf")

        if os.path.exists("cnvkit.pdf"):
            toc.append([1, "Copy Number Plots", len(doc) + offset])
            doc.insert_file("cnvkit.pdf")

        toc.append([1, "Mutational signatures", len(doc) + offset])
        if len(self.files["signatures"]) == 1:
            doc.insert_file(self.files["signatures"][0])
        else:
            for index, s in enumerate(self.files["signatures"]):
                toc.append([2, f"Sample {index}", len(doc) + offset])
                doc.insert_file(s)

        if self.show_therapies and os.path.exists("reference_list.pdf"):
            toc.append([1, "References", len(doc) + offset])
            doc.insert_file("reference_list.pdf")

        self.build_toc(toc)
        doc.insert_file("toc.pdf", start_at=1)

        total_pages: int = len(doc)
        for i, p in enumerate(doc.pages()):
            if i > 0:
                self._add_page_numbers(p, total_pages)
        doc.set_toc(toc)
        if not Path(self.filename).is_absolute():
            doc.save(f"{self.original_dir}/{self.filename}")
        else:
            doc.save(self.filename)
