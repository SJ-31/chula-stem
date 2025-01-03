#!/usr/bin/env ipython
import datetime
import os

import polars as pl
import pymupdf
from pymupdf import Document, Page, mupdf
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import cm, inch
from reportlab.platypus import (
    Paragraph,
    SimpleDocTemplate,
    Table,
)

import chula_stem.report.format as fr
from chula_stem.databases import get_therapy_df
from chula_stem.plotting import plot_cnvkit
from chula_stem.report import ResultsReport, add_pstyles
from chula_stem.report.spec import (
    BOLD_FONT,
    FONT,
    STYLE,
    cnv_style,
    reference_list_style,
    repeat_style,
    signature_style,
    snp_style,
    therapy_style,
)
from chula_stem.report.utils import style_cells


class VariantCallingReport(ResultsReport):
    def __init__(
        self,
        filename: str,
        metadata: dict,
        pandrugs2_cache: str = "",
        civic_cache: str = "",
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
    ) -> None:
        super().__init__(
            filename=filename,
            tmpdir=tmpdir,
            font=font,
            bold_font=bold_font,
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
        self.metadata: dict = metadata
        self.civic_cache: str = civic_cache
        self.header_y: int = cm
        self.pandrugs2_cache: str = pandrugs2_cache
        self.data: dict = {"relevant": {}, "nonrelevant": {}, "all": {}}
        self.files: dict = {"signatures": []}
        self.table_styles: dict = {
            "small": snp_style(),
            "cnv": cnv_style(),
            "repeat": repeat_style(),
        }
        self.table_styles["sv"] = self.table_styles["small"]
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
                    "Tandem Repeats",
                    "Copy Number Variants",
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
                spec["gene_col"] = "All Genes"
                spec["separator"] = ","
            elif type == "repeat":
                spec["gene_col"] = "Affected Gene"
            variant_spec.append(spec)

        self.data["therapies"], self.data["study_references"] = fr.therapy_fmt(
            therapy_df, variant_spec
        )
        self.stats["counts"]["Available Therapies"] = self.data["therapies"].shape[0]
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
        # TODO: front page will be a table of two columns, one for patient info, the other sample info
        fontsize = 12
        widths = [100, 150]
        styles = {
            0: ParagraphStyle("keys", fontName=BOLD_FONT, fontSize=fontsize),
            None: ParagraphStyle("values", fontName=FONT, fontSize=fontsize),
        }
        patient_data = [
            ["Name", self.metadata.get("patient_name", "-")],
            ["Gender", self.metadata.get("gender", "-")],
            ["Date of Birth", self.metadata.get("patient_birthdate", "-")],
            ["Physician", self.metadata.get("physician", "-")],
            ["Diagnosis", self.metadata.get("diagnosis", "-")],
            ["Hospital", self.metadata.get("hospital", "-")],
        ]
        sample_data = [
            ["ID", self.metadata.get("id", "-")],
            ["Collection Date", self.metadata.get("collection_date", "-")],
            ["Type", self.metadata.get("sample_type", "-")],
        ]
        patient_data = add_pstyles(patient_data, styles, nested=True)
        sample_data = add_pstyles(sample_data, styles, nested=True)
        patient_table = Table(patient_data, colWidths=widths)
        sample_table = Table(sample_data, colWidths=widths)
        header = ["Patient Information", "", "Sample Information"]
        header = add_pstyles(
            header,
            ParagraphStyle("header", fontName=BOLD_FONT, fontSize=1.2 * fontsize),
        )
        header_styles = style_cells(
            (0, 0), ncols=3, nrows=1, align="CENTER", background=colors.whitesmoke
        )
        cell_styles = style_cells((0, 1), background=colors.whitesmoke, valign="TOP")
        padding = ["", "", ""]
        full = Table(
            [header, padding, [patient_table, "", sample_table]],
            colWidths=[200, 100, 200],
            style=cell_styles + header_styles,
        )
        doc = SimpleDocTemplate(
            "front_page.pdf",
            pagesize=A4,
            rightMargin=self.spec.get("right_margin", inch),
            leftMargin=self.spec.get("left_margin", inch),
            topMargin=self.spec.get("top_margin", inch),
            bmargin=self.spec.get("bottom_margin", inch),
        )
        doc.build([full])

    def build_toc(self, toc: list) -> None:
        """Build a toc page

        :param toc: toc list of the form used by pymupdf
        """
        toc_fontsize = 8
        cells = []
        header = add_pstyles(
            ["Table of Contents", "Page"],
            {
                None: ParagraphStyle(
                    "header", fontName=BOLD_FONT, fontSize=1.2 * toc_fontsize
                ),
                2: ParagraphStyle(
                    "header",
                    fontName=BOLD_FONT,
                    fontSize=1.2 * toc_fontsize,
                    alignment=2,
                ),
            },
        )
        for level, text, page in toc:
            text_style = ParagraphStyle(
                "text",
                fontName=FONT,
                fontSize=toc_fontsize,
                firstLineIndent=(level - 1) * 20,
            )
            num_style = ParagraphStyle(
                "num", fontName=BOLD_FONT, fontSize=toc_fontsize, alignment=2
            )
            t = Paragraph(text, text_style)
            n = Paragraph(str(page), num_style)
            cells.append([t, n])
        cells = [header, ["", "", ""]] + cells
        toc_table = Table(cells, colWidths=[150, 50], rowHeights=10)
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
        full = Table([[toc_table, "", stats_table]], colWidths=[200, 50, 200])
        doc.build([full])

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
        # TODO: create universal footer/header spec for tables
        # TODO: create universal styles for variant tables
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
        self.build_table(
            self.spec,
            therapy_style(),
            self.data["therapies"],
            "Relevant therapies",
            "Relevant therapies (Continued)",
            table_decorator,
            "therapies.pdf",
        )
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

    # TODO: add in the front, back pages
    def merge(self) -> None:
        toc: list = []
        doc: Document = pymupdf.open()

        offset = 2  # Account for front page and TOC page
        doc.insert_file("front_page.pdf")
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

        toc.append([1, "References", len(doc) + offset])
        doc.insert_file("reference_list.pdf")

        self.build_toc(toc)
        doc.insert_file("toc.pdf", start_at=1)

        total_pages: int = len(doc)
        for p in doc.pages():
            self._add_page_numbers(p, total_pages)
        doc.set_toc(toc)
        doc.save(self.filename)
