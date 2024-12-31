#!/usr/bin/env ipython

# Formatting functions for the report
# Functions ending in the suffix "_fmt" produce formatted dataframes for the
#    report builder. They must return a tuple of "all", "relevant", "nonrelevant"
#    the latter two being used directly for report tables
#    In addition to cleaning up the output with string operations and renaming cols,
#       these fns will add relevant links

import polars as pl
from chula_stem.report.spec import URL, Rename
from chula_stem.utils import add_loc, empty_string2null, read_facets_rds
import os
import re
from chula_stem.databases import add_therapy_info
import polars.selectors as cs
from chula_stem.callset_qc import IMPACT_MAP

TUMOR_KEYWORDS = ["cancer", "leukemia", "carcinoma", "lymphoma"]


def get_genes(file_spec: list[dict]) -> list:
    unique_genes: set = set()
    for spec in file_spec:
        df = pl.read_csv(
            spec["file"],
            separator=spec.get("separator", "\t"),
            null_values=["NA", "."],
            infer_schema_length=None,
        )
        gene_col: pl.Series = df[spec["column"]].drop_nulls()
        if spec.get("is_list"):
            cur_genes: list = [
                g.strip()
                for genes in gene_col.str.split(spec.get("is_list_separator", ";"))
                for g in genes
            ]
        else:
            cur_genes = list(gene_col)
        unique_genes |= set(cur_genes)
    return list(unique_genes)


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
    relevant = with_cn.filter(pl.col("Known/predicted Dosage-sensitive Genes") != "-")
    nonrelevant = with_cn.filter(
        pl.col("Known/predicted Dosage-sensitive Genes") == "-"
    )
    return with_cn, relevant, nonrelevant


def vep_fmt(
    vep_path: str,
    variant_class: str,
) -> tuple:
    """Format and filter vep output into a dataframe with values ready to write
    into a reportlab table
    """
    var_col: str = "Database Name"  # New column for 'Existing_variation'
    df: pl.DataFrame = pl.read_csv(
        vep_path, separator="\t", null_values=["NA", "."], infer_schema_length=None
    ).with_columns(
        pl.concat_str([pl.col("Loc"), pl.col("Feature")], separator="|").alias("VAR_ID")
    )
    if variant_class == "sv":
        rename: dict = Rename.sv
    else:
        rename: dict = Rename.snp
    wanted_cols: list = list(rename.values())
    df = (
        df.drop("Gene")
        .rename(rename)
        .with_columns(
            pl.col("HGVS").str.extract(r".*:(.*)$", 1),
            pl.col("Variant Type").str.replace("_variant$", ""),
            (pl.col(rename["VAF"]).cast(pl.Float64) * 100).round(3),
        )
        .with_columns(
            cs.by_dtype(pl.String)
            .fill_null("-")
            .str.replace_all("&", ", ", literal=True)
            .str.replace_all("_", " ", literal=True),
        )
        .filter(pl.col("Gene") != "-")
        .group_by(pl.col(["Gene", "Locus"]))
        .agg(pl.all().first())
    confident = (
        df.filter((pl.col(var_col).is_not_null()) & (pl.col("SOMATIC") == "1"))
        .with_columns(impact_score=pl.col("IMPACT").replace_strict(IMPACT_MAP))
        .sort(pl.col("impact_score"), descending=True)
    )
    others = df.filter(~pl.col("VAR_ID").is_in(confident["VAR_ID"]))
    relevant = confident.select(wanted_cols)
    nonrelevant = others.select(wanted_cols)
    return df, relevant, nonrelevant


def format_therapy_sources(sources: str, source_lookup: dict) -> list[str]:
    def helper(source):
        if find := re.findall(r"\[(.*)\]\((.*)\)", source):
            match = find[0]
            with_link = add_link(match[0], match[1], underline="yes")
        else:
            return source

        num = source_lookup.get(with_link, f"[{len(source_lookup) + 1}]")
        source_lookup[with_link] = num
        return num

    return [helper(s) for s in sources if s and s != "NA"]


def therapy_fmt(
    therapy_df: pl.DataFrame,
    variant_spec: list[dict],
) -> tuple[pl.DataFrame, pl.DataFrame]:
    def db_link_fn(name, link) -> str:
        params = {"underline": "yes", "underlinecolor": "#5e81ac"}
        if name == "Civic":
            id = link.replace("https://civicdb.org/evidence/", "")
            return add_link(f"Civic:{id}", link, **params)
        else:
            return add_link(name, link, **params)

    source_dict: dict = {}
    dfs: list = []
    for v in variant_spec:
        gene_col: str = v["gene_col"]
        current: pl.DataFrame = v["df"].select(gene_col)
        if v.get("is_list"):
            current = (
                current.with_columns(pl.col(gene_col).str.split(v["separator"]))
                .explode(gene_col)
                .with_columns(pl.col(gene_col).str.strip_chars())
            )
        type: str = v["type"]
        genes: set = set(current[gene_col])
        df = therapy_df.filter(pl.col("gene").is_in(genes)).with_columns(
            pl.lit(type).alias("type")
        )
        dfs.append(df)
    wanted_cols = [
        "gene",
        "disease",
        "type",
        "source",
        "therapies",
        "db",
        "db_link",
    ]

    therapy_df = (
        pl.concat(dfs)
        .select(wanted_cols)
        .explode("therapies")
        .explode("disease")
        .unique()
        .with_columns(
            pl.col("disease").str.to_lowercase().str.replace_all("_", " "),
        )
        .filter(  # Make sure that the therapies reported are relevant to cancer
            (pl.col("db") == "PanDrugs2")
            | pl.col("disease").str.contains_any(TUMOR_KEYWORDS)
        )
        .pipe(empty_string2null)
        .with_columns(
            pl.col("disease")
            .str.replace_many({"cancer": "", "clinical": ""})
            .map_elements(str.capitalize, return_dtype=pl.String),
            pl.struct(["db", "db_link"])
            .map_elements(
                lambda x: db_link_fn(x["db"], x["db_link"]),
                return_dtype=pl.String,
            )
            .alias("db_link"),
        )
        .group_by("therapies")
        .agg(pl.all().unique())
        .with_columns(
            pl.col("source").map_elements(
                lambda x: format_therapy_sources(x, source_dict),
                return_dtype=pl.List(pl.String),
            ),
        )
        .with_columns(
            pl.col("db_link").list.head(5).list.join(", "),
            pl.col(["gene", "disease", "type", "source", "db"]).list.join(", "),
            pl.col("therapies")
            .str.split(":")
            .list.to_struct(fields=["therapies", "PubChemId"]),
        )
        .unnest("therapies")
        .filter(pl.col("therapies").str.contains("[A-Z-a-z]*"))
        .with_columns(
            pl.col("PubChemId")
            .map_elements(
                lambda x: (
                    add_link(x, f"{URL.pubchem}/{x}", underline="yes")
                    if x != "NA"
                    else x
                ),
                return_dtype=pl.String,
            )
            .str.replace("NA", "-"),
        )
        .rename(Rename.therapy)
        .select(list((Rename.therapy).values()))
        .pipe(empty_string2null)
        .fill_null("-")
    )
    source_df = pl.DataFrame(
        {"Number": source_dict.values(), "Source": source_dict.keys()}
    )
    return therapy_df, source_df


def format_signature_link(sig: str, link: str) -> str:
    return add_link(sig, link, underline="yes")


def make_signature_table(slice: pl.DataFrame, ref: pl.DataFrame) -> pl.DataFrame:
    """Format a row of the signature df (a single sample) into a df to present with
    reportlab

    The signature df has four columns: Samples, Signatures, Frequency and m
    Signatures and Frequency are list columns of the same length, containing the
    names and frequencies of signatures kept in the sample in their respective orders

    :returns: a df of formatted data
    """
    wanted_cols = [
        "Signature",
        "Count",
        "Frequency",
        "Collection",
        "Proposed Aetiology",
    ]
    total: int = slice["m"][0]
    df = (
        (
            slice.explode(["Signatures", "Count", "Frequency"])
            .join(ref, left_on="Signatures", right_on="Signature", how="left")
            .rename(lambda x: x.replace("_", " "))
            .with_columns(
                pl.struct(s="Signatures", l="Link")
                .map_elements(
                    lambda x: format_signature_link(x["s"], x["l"]),
                    return_dtype=pl.String,
                )
                .alias("Signature")
            )
        )
        .select(wanted_cols)
        .rename({"Count": f"Count (Total: {total})"})
    )
    return df


def sigprofiler_fmt(
    soln_activities_path: str,
    cosmic_reference: str,
    abs_threshold: int = 100,
    rel_threshold: float = 0.25,
    excluded_signatures: str = "",
) -> list[pl.DataFrame]:
    """Format SigProfilerAssignment results for reportlab

    :param soln_activities_path: Path to SigProfiler "Assignment_Solution_Activities.txt" file
    :param cosmic_reference: Path to csv file with COSMIC signature metadata

        This function performs some basic filtering on signatures
    :param abs_threshold: Minimum absolute number of mutations required to be considered
    :param rel_threshold: Minimum percentage of mutational activity to be considered
    :param excluded_signatures: Path to file containing signatures to ignore
        Signatures in this file are spearated by newlines

    :returns: a list of (n mutations in sample, sample dataframe)
        In case only one sample was provided to SigProfiler, just take the first element
    """
    if excluded_signatures:
        with open(excluded_signatures, "r") as f:
            excluded: list = f.readlines()
    else:
        excluded: list = []
    df = pl.read_csv(soln_activities_path, separator="\t").select(
        ~cs.by_name(*excluded)
    )
    signatures: list = df.columns[1:]
    sig_cols = pl.col(signatures)
    samples: pl.Series = df["Samples"]
    replace_expr = [
        pl.col(s).replace_strict({True: s, False: None}) for s in signatures
    ]

    m: pl.Series = df.select(sig_cols).sum_horizontal()
    # Total number of signature mutations per sample
    frequencies = (
        df.with_columns(sig_cols / m)
        .unpivot(
            on=signatures,
            index="Samples",
            variable_name="Signatures",
            value_name="Frequency",
        )
        .with_columns(pl.col("Frequency").round(2))
    )
    counts = df.with_columns(sig_cols).unpivot(
        on=signatures, index="Samples", variable_name="Signatures", value_name="Count"
    )

    sums: pl.DataFrame = pl.DataFrame({"Samples": samples, "m": m})
    filtered = (
        (
            df.with_columns(
                (sig_cols > abs_threshold) | ((sig_cols / m) > rel_threshold)
            )
            .with_columns(replace_expr)
            .with_columns(Signatures=pl.concat_list(signatures).list.drop_nulls())
            .select(["Samples", "Signatures"])
            .join(sums, on="Samples")
            .explode("Signatures")
        )
        .join(frequencies, on=["Samples", "Signatures"])
        .join(counts, on=["Samples", "Signatures"])
        .group_by("Samples")
        .agg(pl.col(["Signatures", "Frequency", "Count"]), pl.col("m").first())
    )
    cosmic_data: pl.DataFrame = pl.read_csv(cosmic_reference)
    formatted = [make_signature_table(s, cosmic_data) for s in filtered.iter_slices(1)]
    return formatted
