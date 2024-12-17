import chula_stem.report as rp
import polars as pl
from chula_stem.callset_qc import IMPACT_MAP
from chula_stem.report import Civic, PanDrugs2, TherapyDB

small_path = (
    "/home/shannc/Bio_SDD/chula-stem/tests/vep/7-patient_10-VEP_small_1.tsv"
)
df: pl.DataFrame = pl.read_csv(
    small_path, separator="\t", null_values=["NA", "."], infer_schema_length=None
).with_columns(
    pl.concat_str([pl.col("Loc"), pl.col("Feature")], separator="|").alias(
        "VAR_ID"
    )
)

##
cache: str = "/home/shannc/Bio_SDD/chula-stem/tests/civic.json"
pd_cache: str = "/home/shannc/Bio_SDD/chula-stem/tests/pandrugs2.json"
therapy_dbs: list[TherapyDB] = [Civic(cache), PanDrugs2(pd_cache)]
gene_list = list(df["SYMBOL"].unique())

gene_list: list = gene_list[:10]  # Truncate for testing purposes
gene_list: list = [
    "BRAF",
    "CEBPA",
    "CDK6",
    "HMGA2",
    "CCNE1",
    "FGFR1",
    "BRCA1",
    "MGMT",
    "ALK",
    "EGFR",
    "KRAS",
]

## * Get therapy data
temp: list[pl.DataFrame] = []
for db in therapy_dbs:
    find_info: pl.DataFrame = db.get_genes(gene_list, True)
    if not find_info.is_empty():
        found = find_info["gene"]
        gene_list = list(filter(lambda x: x not in found, gene_list))
        temp.append(find_info.select(TherapyDB.shared_cols))
    if not gene_list:
        break
all_drug_info: pl.DataFrame = pl.concat(temp)
with_therapeutics = df.join(
    all_drug_info, how="left", left_on="SYMBOL", right_on="gene"
).unique("SYMBOL")
confident = (
    with_therapeutics.filter(
        (pl.col("Existing_variation").is_not_null()) & (pl.col("SOMATIC") == "1")
    )
    .with_columns(impact_score=pl.col("IMPACT").replace_strict(IMPACT_MAP))
    .sort(pl.col("impact_score"), descending=True)
)
others = confident.filter(~pl.col("VAR_ID").is_in(confident["VAR_ID"]))

# The table for renaming VEP columns for the report
VTABLE_RENAME: dict = {
    "SYMBOL": "Gene",
    "VAF": "Variant Allele Frequency",
    "Alt_depth": "Variant Read Support",
    "Loc": "Locus",
    "HGVSc": "HGVS",
    "Existing_variation": "Database Name",
    "Consequence": "Variant Type",
    "CLIN_SIG": "ClinVar",
}
## * Report generation
VTABLE_COL_WIDTHS: dict = {
    "Gene": 70,
    "Variant Allele Frequency": 60,
    "Variant Read Support": 50,
    "Locus": 80,
    "HGVS": 90,
    "Variant Type": 90,
    "ClinVar": 100,
}

from reportlab.lib import colors, units
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfgen import canvas
from reportlab.platypus import LongTable, Paragraph, SimpleDocTemplate, Spacer

# (0, 0) is the bottom left of the page by default


# It'll be best to wrap this in a class
class ResultsReport:
    def __init__(self, filename: str) -> None:
        self.canv: canvas.Canvas = canvas.Canvas(filename)
        pass

    # def front_page(self) -> None:

    def draw_variant_table(
        self, variants: pl.DataFrame, width: int = 900
    ) -> None:
        res = (
            variants.drop("Gene")
            .rename(VTABLE_RENAME)
            .select(VTABLE_RENAME.values())
        )
        to_data = [row for row in res.iter_rows()]


# Flowables
# - Flowable.drawOn(canvas, x, y)

# Documents https://docs.reportlab.com/reportlab/userguide/ch5_platypus/#documents-and-templates
# Create a list of flowables and pass it to Doc.build()
import polars.selectors as cs

res = (
    with_therapeutics.drop("Gene")
    .rename(VTABLE_RENAME)
    .select(VTABLE_RENAME.values())
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

# Using Paragraphs to format text for word wrap
column_style: ParagraphStyle = ParagraphStyle("cols")

to_data = rp.add_pstyles(res, column_style)
cols = rp.add_pstyles(res.columns, column_style)
to_data.insert(0, cols)

# repeatRows=1 repeats the first row at every split
T = LongTable(to_data, colWidths=list(VTABLE_COL_WIDTHS.values()), repeatRows=1)
ncols: int = len(res.columns)

# Dummy styles for now
col_styles = rp.style_cells(
    (0, 0),
    ncols,
    1,
    textcolor=colors.red,
    fontsize=13,
    underline=(3, colors.black),
    background=colors.lightgrey,
)
style_all = rp.style_cells((0, 1), background=colors.lightcyan, valign="TOP")

# Use a doc template to enable the table to be split across multiple pages
pdf = SimpleDocTemplate(
    "/home/shannc/Bio_SDD/chula-stem/test.pdf", rightMargin=3, leftMargin=3
)
T.setStyle(col_styles + style_all)
loc: tuple = (10, 10)
# splits = T.splitOn(pdf, 500, 100)
elements = []
elements.append(T)
pdf.build(elements)
