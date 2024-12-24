#!/usr/bin/env ipython

# Formatting functions for the report
# Functions ending in the suffix "_fmt" produce formatted dataframes for the
#    report builder. They must return a tuple of "all", "relevant", "nonrelevant"
#    the latter two being used directly for report tables
#    In addition to cleaning up the output with string operations and renaming cols,
#       these fns will add relevant links

import polars as pl
from chula_stem.report.spec import URL, Rename
from chula_stem.utils import add_loc, read_facets_rds
import os
from chula_stem.databases import add_therapy_info
import polars.selectors as cs
from chula_stem.callset_qc import IMPACT_MAP

TUMOR_KEYWORDS = ["cancer", "leukemia", "carcinoma", "lymphoma"]


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


def dbvar_link(x):
    link = add_link(x, f"{URL.dbvar}/{x}", underline="yes")
    return f"dbVar:{link}"


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


def classify_cnv(
    classifycnv_path,
    facets_path: str = "",
    cnvkit_path: str = "",
):
    cnv = pl.read_csv(
        classifycnv_path, separator="\t", null_values="NA", infer_schema_length=None
    ).with_columns(
        pl.struct(s="source", acc="accession")
        .map_elements(
            lambda x: (dbvar_link(x["acc"]) if x["acc"] != "NA" else x["s"]),
            return_dtype=pl.String,
        )
        .alias("source")
    )
    wanted_cols = ["Locus"] + list(Rename.cnv.values())
    cnv = (
        cnv.with_columns(
            pl.concat_str(["Start", "End"], separator="-").alias("range"),
        )
        .with_columns(
            pl.concat_str(["Chromosome", "range"], separator=":").alias("Locus"),
        )
        .rename(Rename.cnv)
    ).select(wanted_cols)
    cn_col = "Estimated Copy Number"
    wanted_cols.insert(1, cn_col)
    with_cn = (
        get_copy_number(cnv, facets_path, cnvkit_path, cn_col_name="cn")
        .rename({"cn": cn_col})
        .select(wanted_cols)
    )
    relevant = with_cn.filter(
        pl.col("Known/predicted dosage-sensitive genes").is_not_null()
    )
    nonrelevant = with_cn.filter(
        pl.col("Known/predicted dosage-sensitive genes").is_null()
    )
    return with_cn, relevant, nonrelevant


def vep_fmt(
    vep_path: str,
    tmpdir: str,
    variant_class: str,
    civic_cache: str,
    pandrugs2_cache: str,
) -> tuple:
    """Format and filter vep output into a dataframe with values ready to write
    into a reportlab table
    """
    r_out = f"{tmpdir}/{variant_class}_relevant.parquet"
    nr_out = f"{tmpdir}/{variant_class}_non-relevant.parquet"
    all_out = f"{tmpdir}/{variant_class}_all.parquet"
    if os.path.exists(r_out) and os.path.exists(nr_out) and os.path.exists(all_out):
        return tuple([pl.read_parquet(p) for p in [all_out, r_out, nr_out]])
    var_col: str = "Database Name"  # New column for 'Existing_variation'
    with_therapeutics: pl.DataFrame = add_therapy_info(
        vep_path, civic_cache, pandrugs2_cache
    )
    wanted_cols: list = list(Rename.sv_snp.values())
    with_therapeutics = (
        with_therapeutics.drop("Gene")
        .rename(Rename.sv_snp)
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
    all = with_therapeutics
    relevant = confident.select(wanted_cols)
    nonrelevant = others.select(wanted_cols)
    for df, path in zip([all, relevant, nonrelevant], [r_out, nr_out, all_out]):
        df.write_parquet(path)
    return all, relevant, nonrelevant
