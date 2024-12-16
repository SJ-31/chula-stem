import chula_stem.report as rp
import polars as pl
from chula_stem.callset_qc import IMPACT_MAP
from chula_stem.report import Civic, PanDrugs2, TherapyDB

small_path = "/home/shannc/Bio_SDD/chula-stem/tests/vep/small.tsv"
df: pl.DataFrame = pl.read_csv(
    small_path, separator="\t", null_values="NA", infer_schema_length=None
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
)
confident = (
    with_therapeutics.filter(
        (pl.col("Existing_variation").is_not_null()) & (pl.col("SOMATIC") == "1")
    )
    .with_columns(impact_score=pl.col("IMPACT").replace_strict(IMPACT_MAP))
    .sort(pl.col("impact_score"), descending=True)
)
others = confident.filter(~pl.col("VAR_ID").is_in(confident["VAR_ID"]))

vtable_rename: dict = {
    "SYMBOL": "Gene",
    "VAF": "Variant Allele Frequency",
    "Alt_depth": "Variant Read Support",
    "Loc": "Locus",
    "HGVSc": "HGVS",
    "Existing_variation": "Database Name",
    "Consequence": "Variant Effect",
    "CLIN_SIG": "ClinVar",
}
