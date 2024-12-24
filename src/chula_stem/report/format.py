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
import re
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


def get_clingen_link(clingen_str) -> str:
    def get_sym(clingen_link: str) -> str:
        find = re.findall(r".*sym=(.*)\&.*", clingen_link)
        if not find:
            return ""
        return find[0]

    symbol = get_sym(clingen_str)
    if not symbol:
        return "NA"
    return add_link(symbol, clingen_str, underline="yes")


def msisensor_pro_fmt(msisensor_path: str) -> tuple:
    wanted_cols = ["Locus"] + list(Rename.repeat.values())
    df = pl.read_csv(
        msisensor_path, separator="\t", null_values="NA", infer_schema_length=None
    )
    df = (
        (
            add_loc(df, start_col="Start", end_col="End")
            .with_columns(
                pl.col("ClinGen_report").map_elements(
                    get_clingen_link, return_dtype=pl.String
                )
            )
            .with_columns(
                pl.struct(s="source", acc="accession")
                .map_elements(
                    lambda x: (
                        dbvar_link(x["acc"])
                        if x["acc"] != "NA" and x["acc"]
                        else x["s"]
                    ),
                    return_dtype=pl.String,
                )
                .alias("source")
            )
        )
        .rename(Rename.repeat)
        .fill_null("-")
    )
    relevant = df.filter(pl.col("ClinGen") != "-").select(wanted_cols)
    nonrelevant = df.filter(pl.col("ClinGen") != "-").select(wanted_cols)
    return df, relevant, nonrelevant


def classify_cnv_fmt(
    classifycnv_path,
    facets_path: str = "",
    cnvkit_path: str = "",
) -> tuple:
    """
    Format cnv data from ClassifyCNV for the report. Includes merging with cnvkit
    and facets to determine show estimated copy number
    """
    cnv = pl.read_csv(
        classifycnv_path, separator="\t", null_values="NA", infer_schema_length=None
    ).with_columns(
        pl.struct(s="source", acc="accession")
        .map_elements(
            lambda x: (
                dbvar_link(x["acc"]) if x["acc"] != "NA" and x["acc"] else x["s"]
            ),
            return_dtype=pl.String,
        )
        .alias("source"),
        pl.col("ClinGen_report").map_elements(get_clingen_link, return_dtype=pl.String),
    )
    # "group_by" operation required because data was exploded to become longer on
    #   dosage-sensitive genes in the R script
    wanted_cols = ["Locus"] + list(Rename.cnv.values())
    cnv = (
        add_loc(cnv, start_col="Start", end_col="End", chr_col="Chromosome")
        .rename(Rename.cnv)
        .group_by("Locus")
        .agg((~cs.by_name("ClinGen")).unique().first(), pl.col("ClinGen"))
        .with_columns(pl.col("ClinGen").list.join(", "))
    )
    cn_col = "Estimated Copy Number"
    wanted_cols.insert(1, cn_col)
    with_cn = (
        get_copy_number(cnv, facets_path, cnvkit_path, cn_col_name="cn")
        .rename({"cn": cn_col})
        .select(wanted_cols)
    ).fill_null("-")
    relevant = with_cn.filter(
        pl.col("Known/predicted Dosage-sensitive Genes") != "-"
    )
    nonrelevant = with_cn.filter(
        pl.col("Known/predicted Dosage-sensitive Genes") == "-"
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
    for df, path in zip([relevant, nonrelevant, all], [r_out, nr_out, all_out]):
        df.write_parquet(path)
    return all, relevant, nonrelevant


def therapy_fmt(parquet_path):
    # TODO: must include signatures
    db_link_fn = lambda x, y: add_link(x, y, underline="yes", underlinecolor="#5e81ac")
    df = (
        (
            (
                pl.read_parquet(parquet_path)
                .select(
                    pl.col(["Gene", "disease", "source", "therapies", "db", "db_link"])
                )
                .explode("disease")
            )
            .with_columns(
                pl.col("disease").str.to_lowercase().str.replace_all("_", " "),
                pl.struct(["db", "db_link"])
                .map_elements(
                    lambda x: db_link_fn(x["db"], x["db_link"]),
                    return_dtype=pl.String,
                )
                .alias("db_link"),
            )
            .filter(  # Make sure that the therapies reported are relevant to cancer
                (pl.col("db") == "pandrugs2")
                | pl.col("disease").str.contains_any(TUMOR_KEYWORDS)
            )
            .with_columns(
                pl.col("disease")
                .str.replace_many({"cancer": "", "clinical": ""})
                .str.replace("", "unspecified")
                .map_elements(str.capitalize, return_dtype=pl.String)
            )
        )
        .explode("therapies")
        .with_columns(
            pl.col("therapies")
            .str.split(":")
            .list.to_struct(fields=["therapies", "PubChemId"])
        )
        .unnest("therapies")
        .with_columns(
            pl.col("PubChemId").map_elements(
                lambda x: (
                    add_link(x, f"{URL.pubchem}/{x}", underline="yes")
                    if x != "NA"
                    else x
                ),
                return_dtype=pl.String,
            ),
            pl.lit("Gene variation").alias("type"),
        )
        .filter(pl.col("therapies").str.contains("[A-Z-a-z]*"))
        .rename(Rename.therapy)
    ).select(list((Rename.therapy).values()))
    relevant = df.filter(
        (pl.col("Relevant cancers") != "Unspecified") & (pl.col("Study source") != "NA")
    )
    nonrelevant = df.filter(~pl.col("Therapy").is_in(relevant["Therapy"]))
    return df, relevant, nonrelevant
