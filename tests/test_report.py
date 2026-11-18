from pathlib import Path

import pytest
from pyhere import here

# Documents https://docs.reportlab.com/reportlab/userguide/ch5_platypus/#documents-and-templates
# Create a list of flowables and pass it to Doc.build()
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import BaseDocTemplate
from sarek_plus import VariantCallingReport, entry_point
from sarek_plus.report.format import (
    AscatResult,
    CNVKitResult,
    CopyNumberResult,
    FreecResult,
)

STYLES = getSampleStyleSheet()
PWIDTH = A4[0]
PHEIGHT = A4[1]

TEST_DIR = here("tests")


# (0, 0) is the bottom left of the page by default

# Flowables
# - Flowable.drawOn(canvas, x, y)

# You can work in cm by converting doc.width (or height) with the `units` module like so:
# doc.width/units.cm
# Use a doc template to enable the table to be split across multiple pages
out = "/home/shannc/Bio_SDD/chula-stem/test.pdf"
# For doc coordinates, (0, 0) is the bottom left corner


def default_shapes(canvas, doc: BaseDocTemplate):
    canvas.saveState()
    bottom_line_h = inch * 1
    top_line_h = PHEIGHT - inch * 1
    canvas.line(0, bottom_line_h, PWIDTH, bottom_line_h)
    canvas.line(0, top_line_h, PWIDTH, top_line_h)
    canvas.restoreState()


@pytest.mark.skip(reason="Done")
def test_full(front_page_only=False):
    from sarek_plus.report.variant_calling_report import VariantCallingReport

    meta = {
        "patient_name": "John Doe",
        "gender": "Male",
        "patient_birthdate": "1990-05-14",
        "physician": "Dr. Alice Smith",
        "diagnosis": "Hypertension",
        "hospital": "Green Valley Medical Center",
        "id": "SMP-123456",
        "collection_date": "2024-12-01",
        "sample_type": "Blood",
    }
    R = VariantCallingReport(
        "/home/shannc/Bio_SDD/chula-stem/tests/report_full.pdf",
        metadata=meta,
        vep_small=small_path,
        vep_sv=sv_pat,
        classify_cnv="/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-classify_cnv-CR.tsv",
        msisensor_pro="/home/shannc/Bio_SDD/chula-stem/tests/msisensor/4-null-CR.tsv",
        sigprofiler="/home/shannc/Bio_SDD/chula-stem/tests/6-null-SigProfilerAssignment/Activities/Assignment_Solution_Activities.txt",
        cosmic_signatures="/home/shannc/Bio_SDD/chula-stem/nextflow/config/cosmic_signatures_v3.4-2024-12-26.csv",
        tmpdir="/home/shannc/Bio_SDD/chula-stem/tests/report_tmp",
    )
    if front_page_only:
        R.build_front_page()
    else:
        R.build()
        R.merge()


@pytest.mark.skip(reason="Done")
def test_front_page():
    test_full(True)


spec: list = [
    {
        "file": "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-classify_cnv-CR.tsv",
        "is_list": True,
        "is_list_separator": ",",
        "column": "All protein coding genes",
    },
    {
        "file": "/home/shannc/Bio_SDD/chula-stem/tests/msisensor/5-P1-Msisensor_unstable.tsv",
        "column": "gene_name",
    },
    {
        "file": "/home/shannc/Bio_SDD/chula-stem/tests/vep/sv.tsv",
        "column": "SYMBOL",
    },
]


@pytest.mark.skip(reason="Done")
def test_read_json():
    path = "/home/shannc/Bio_SDD/chula-stem/tests/report_params.json"
    entry_point.callback(
        report_type="variant_calling",
        specification=path,
        tmpdir="/home/shannc/Bio_SDD/chula-stem/tests/report_tmp",
        output="/home/shannc/Bio_SDD/chula-stem/report_full.pdf",
    )


def test_cn_intervals():
    from intervaltree import IntervalTree
    from sarek_plus.report.format import IntervalAnnotation

    list_A = [
        (1487442, 1646733, IntervalAnnotation({"A": 3})),
        (2315891, 2390603, IntervalAnnotation({"A": 0})),
        (3335548, 3729032, IntervalAnnotation({"A": 2})),
        (4450839, 4693528, IntervalAnnotation({"A": 3})),
        (4842924, 4964773, IntervalAnnotation({"A": 2})),
        (5426223, 5813473, IntervalAnnotation({"A": 3})),
        (5921762, 6250450, IntervalAnnotation({"A": 1})),
        (6728865, 7058741, IntervalAnnotation({"A": 1})),
        (7125036, 7529787, IntervalAnnotation({"A": 1})),
        (8128563, 8575971, IntervalAnnotation({"A": 2})),
        (9172157, 9489875, IntervalAnnotation({"A": 1})),
    ]

    list_B = [
        (1224145, 1615131, IntervalAnnotation({"B": 3})),
        (2843649, 2914206, IntervalAnnotation({"B": 2})),
        (2966237, 3274514, IntervalAnnotation({"B": 1})),
        (4270616, 4375535, IntervalAnnotation({"B": 0})),
        (5201341, 5380524, IntervalAnnotation({"B": 3})),
        (6810584, 7122450, IntervalAnnotation({"B": 1})),
        (7326273, 7671241, IntervalAnnotation({"B": 3})),
        (7872368, 7941357, IntervalAnnotation({"B": 3})),
        (8553566, 8762971, IntervalAnnotation({"B": 1})),
        (9011874, 9156036, IntervalAnnotation({"B": 1})),
    ]

    list_C = [
        (1723601, 1873570, IntervalAnnotation({"C": 4})),
        (1877102, 2006003, IntervalAnnotation({"C": 2})),
        (2182413, 2233950, IntervalAnnotation({"C": 2})),
        (2581059, 2641908, IntervalAnnotation({"C": 0})),
        (2824929, 3198538, IntervalAnnotation({"C": 1})),
        (3799468, 3884655, IntervalAnnotation({"C": 4})),
        (4489631, 4827470, IntervalAnnotation({"C": 4})),
        (4903221, 5351476, IntervalAnnotation({"C": 2})),
        (5871323, 6094079, IntervalAnnotation({"C": 2})),
        (6374150, 6461606, IntervalAnnotation({"C": 2})),
        (7321509, 7424066, IntervalAnnotation({"C": 4})),
        (7780356, 7913951, IntervalAnnotation({"C": 3})),
        (8018115, 8178308, IntervalAnnotation({"C": 1})),
        (8210497, 8328830, IntervalAnnotation({"C": 2})),
        (8717225, 8965330, IntervalAnnotation({"C": 3})),
    ]

    tree_a = IntervalTree.from_tuples(list_A)
    tree_b = IntervalTree.from_tuples(list_B)
    tree_c = IntervalTree.from_tuples(list_C)

    merged: IntervalTree = tree_a | tree_b | tree_c
    merged.merge_overlaps(data_reducer=lambda x, y: x + y)
    agg, dis = CopyNumberResult.intervals_to_df(merged)


def test_combine_cnv():
    data_dir = TEST_DIR / "data"
    ascat = (
        data_dir
        / "ascat/PHcase_1-tumor_vs_PHcase_1-normal/PHcase_1-tumor_vs_PHcase_1-normal.cnvs.txt"
    )
    cnvkit = (
        data_dir / "cnvkit/PHcase_1-tumor_vs_PHcase_1-normal/PHcase_1-tumor.call.cns"
    )
    freec = (
        data_dir
        / "freec/PHcase_1-tumor_vs_PHcase_1-normal/PHcase_1-tumor_vs_PHcase_1-normal.p.value.txt"
    )
    dfs = [
        AscatResult().read_segmental(ascat),
        FreecResult().read_segmental(freec),
        CNVKitResult().read_segmental(cnvkit),
    ]
    c, d = CopyNumberResult.combine(dfs, add_genes=True)
    print(c)
    print(d)


def gen_report_paths(sample) -> dict[str, Path]:
    sarek_plus = here("output", "nf-work", "sarek_plus_test", "out")
    sarek = here("output", "sarek_pdac_test")
    paths: dict[str, Path] = {
        "optitype": sarek_plus / "optitype" / sample / f"{sample}-optitype_result.tsv",
        "therapy_data_cache": sarek_plus / "gene_data.db",
        "vep_small_somatic": sarek_plus
        / "annotation"
        / sample
        / f"{sample}-tumor_vs_{sample}-normal.tsv",
        "vcf_format_germline": sarek_plus
        / "format_fields"
        / sample
        / f"{sample}-format_fields_germline.tsv",
        "vcf_format_somatic": sarek_plus
        / "format_fields"
        / sample
        / f"{sample}-format_fields_paired.tsv",
        "vep_small_germline": sarek_plus
        / "annotation"
        / sample
        / f"{sample}-normal.tsv",
        "vep_sv_germline": None,
        "vep_sv_somatic": sarek_plus
        / "annotation"
        / sample
        / f"{sample}-sv_concat.tsv",
        "cosmic_signatures": "TODO",
        "classify_cnv": "TODO",
        "sigprofiler": "TODO",
    }
    return {k: v for k, v in paths.items() if v.exists()}


def test_report():
    sample = "PHcase_1"
    sarek_plus = here("output", "nf-work", "sarek_plus_test", "out")
    sarek = here("output", "sarek_pdac_test")
    outdir = here("tests", "sarek_plus_out")
    outfile = outdir / "report.pdf"
    paths = gen_report_paths(sample)
    report = VariantCallingReport(
        filename=outfile,
        other_text={
            "title": "Test WES",
            "disclaimer": "DISCLAIMER",
            "front_page_note": "FRONT_PAGE_NOTE",
            "details": "DETAILS",
        },
        paths=paths,
        metadata={"age": 42, "sex": "YY"},
        tmpdir=outdir / "report_tmp",
        therapy_databases=["open_targets", "civic"],
        download_therapy_data=False,
    )
    report.build()
    report.merge()
