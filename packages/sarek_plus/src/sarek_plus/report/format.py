#!/usr/bin/env ipython

from __future__ import annotations

# Formatting functions for the report
# Functions ending in the suffix "_fmt" produce formatted dataframes for the
#    report builder. They must return a tuple of "all", "relevant", "nonrelevant"
#    the latter two being used directly for report tables
#    In addition to cleaning up the output with string operations and renaming cols,
#       these fns will add relevant links
import re
from collections.abc import Iterable
from copy import deepcopy
from pathlib import Path
from typing import ClassVar, Literal

import numpy as np
import pandera.polars as pa
import polars as pl
import polars.selectors as cs
from attrs import define, field
from chula_stem.utils import (
    Ensembl116,
    add_loc,
    empty_string2null,
    str_split_unique,
)
from intervaltree import Interval, IntervalTree

from sarek_plus.callset_qc import CLINSIG_MAP, IMPACT_MAP, filter_multiallelic
from sarek_plus.databases import TherapyDB
from sarek_plus.mavedb import MaveDB
from sarek_plus.report.spec import URL, Rename

TUMOR_KEYWORDS = ["cancer", "leukemia", "carcinoma", "lymphoma"]


def get_genes(file_spec: list[dict]) -> list:
    unique_genes: set = set()
    for spec in file_spec:
        if spec["file"] is None:
            continue
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
    if "," in x:
        splits = map(lambda y: y.strip(), x.split(","))
        links = [add_link(s, f"{URL.dbvar}/{s}", underline="yes") for s in splits]
        link = ", ".join(links)
    else:
        link = add_link(x, f"{URL.dbvar}/{x}", underline="yes")
    return f"dbVar:{link}"


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


# TODO: look into how to report allele-specific copy number i.e.
# copy number of specific variants


@define
class IntervalAnnotation:
    callers: dict[str, int]

    def __add__(self, anno: IntervalAnnotation) -> IntervalAnnotation:
        return IntervalAnnotation(self.callers | anno.callers)

    def cn(self) -> int | None:
        if self.concordant() and self.callers:
            return next(iter(self.callers.values()))

    def concordant(self) -> bool:
        return len(set(self.callers.values())) == 1


@define
class CopyNumberResult:
    cfg: dict = field(factory=dict)
    segment_schema: ClassVar[pa.DataFrameSchema] = pa.DataFrameSchema(
        {
            "chr": pa.Column(str),
            "start": pa.Column(int),
            "end": pa.Column(int),
            "cn": pa.Column(int),
            "caller": pa.Column(str),
        },
        coerce=True,
    )

    def _read_segmental(self, file: Path) -> pl.LazyFrame:
        raise NotImplementedError()

    def read_segmental(self, file: Path, add_genes: bool = True) -> pl.DataFrame:
        lf: pl.LazyFrame = self._read_segmental(file)
        to_select = list(self.segment_schema.columns.keys())
        if add_genes:
            if "genes" not in lf.collect_schema().names():
                lf = CopyNumberResult.add_gene_data(lf, "chr", "start", "end")
            to_select.extend(["genes", "n_genes"])
        lf = (
            self.segment_schema.validate(lf)
            .with_columns(
                pl.when(pl.col("cn") == 2)
                .then(pl.lit("neutral"))
                .when(pl.col("cn") > 2)
                .then(pl.lit("gain"))
                .otherwise(pl.lit("loss"))
                .alias("cn_category")
            )
            .collect()
        )
        return lf.select(to_select)

    @staticmethod
    def add_gene_data(
        lf: pl.LazyFrame, chr_col: str, start_col: str, end_col: str
    ) -> pl.LazyFrame:
        name_mapping: dict = Ensembl116.name_mapping("ensembl", "gene_name", False)
        gene_coords: pl.LazyFrame = (
            Ensembl116.coords("gene", False)
            .with_columns(pl.col("ensembl").replace(name_mapping))
            .rename({"chr": chr_col, "start": start_col, "end": end_col})
            .lazy()
        )
        lf = (
            lf.with_columns(
                pl.when(pl.col("chr").str.starts_with("chr"))
                .then(pl.col("chr").str.strip_prefix("chr"))
                .otherwise(pl.col("chr"))
            )
            .join_where(
                gene_coords,
                pl.col(chr_col) == pl.col(f"{chr_col}_right"),
                pl.col(start_col) <= pl.col(f"{start_col}_right"),
                pl.col(end_col) >= pl.col(f"{end_col}_right"),
            )
            .drop(cs.ends_with("_right"))
        )
        group_cols = [c for c in lf.collect_schema().names() if c != "ensembl"]
        lf = (
            lf.group_by(group_cols)
            .agg(pl.col("ensembl"))
            .with_columns(pl.col("ensembl").list.unique())
            .rename({"ensembl": "genes"})
            .with_columns(pl.col("genes").list.len().alias("n_genes"))
        )
        return lf

    @staticmethod
    def intervals_to_df(tree: IntervalTree) -> pl.DataFrame:
        tmp_concordant = {"cn": [], "start": [], "end": [], "callers": []}
        tmp_discordant = {"start": [], "end": [], "caller_disagreement": []}
        for el in tree:
            el: Interval
            anno: IntervalAnnotation = el.data
            if anno.concordant():
                tmp_concordant["start"].append(el.begin)
                tmp_concordant["end"].append(el.end)
                tmp_concordant["cn"].append(anno.cn())
                tmp_concordant["callers"].append(",".join(anno.callers))
            else:
                tmp_discordant["start"].append(el.begin)
                tmp_discordant["end"].append(el.end)
                tmp_discordant["caller_disagreement"].append(
                    ",".join([f"{k}={v}" for k, v in anno.callers.items()])
                )

        return pl.DataFrame(tmp_concordant), pl.DataFrame(tmp_discordant)

    @staticmethod
    def combine(
        dfs: list[pl.DataFrame], add_genes: bool = True
    ) -> tuple[pl.DataFrame, pl.DataFrame]:
        """
        Merge copy number results from various callers, returning
        a tuple of concordant, discordant cn regions
        """
        concordant_tmp, discordant_tmp = [], []
        combined: pl.DataFrame = pl.concat(
            [CopyNumberResult.segment_schema.validate(d) for d in dfs],
            how="diagonal_relaxed",
        )
        for chr in combined["chr"].unique():
            cur = combined.filter(pl.col("chr") == chr)
            tree: IntervalTree = IntervalTree.from_tuples(
                [
                    (
                        row["start"],
                        row["end"],
                        IntervalAnnotation({row["caller"]: row["cn"]}),
                    )
                    for row in cur.iter_rows(named=True)
                ]
            )
            tree.merge_overlaps(data_reducer=lambda x, y: x + y)
            conc, disc = CopyNumberResult.intervals_to_df(tree)
            concordant_tmp.append(conc.with_columns(pl.lit(chr).alias("chr")))
            discordant_tmp.append(disc.with_columns(pl.lit(chr).alias("chr")))
        concordant: pl.DataFrame = pl.concat(concordant_tmp, how="diagonal_relaxed")
        discordant: pl.DataFrame = pl.concat(discordant_tmp, how="diagonal_relaxed")
        if add_genes:
            concordant, discordant = [
                CopyNumberResult.add_gene_data(lf.lazy(), "chr", "start", "end")
                .collect()
                .with_columns(pl.col("genes").list.join(", "))
                for lf in (concordant, discordant)
            ]
        return concordant, discordant


@define
class AscatResult(CopyNumberResult):
    schema: pa.DataFrameSchema = pa.DataFrameSchema(
        {
            "chr": pa.Column(pl.String),
            "startpos": pa.Column(pl.Int64),
            "endpos": pa.Column(pl.Int64),
            "nMajor": pa.Column(pl.Int64),
            "nMinor": pa.Column(pl.Int64),
        },
        coerce=True,
    )

    def _read_segmental(self, file: Path) -> pl.LazyFrame:
        "Reads the 'cnvs.txt' file"
        df = (
            pl.scan_csv(file, separator="\t", infer_schema_length=None)
            .pipe(self.schema.validate)
            .rename({"startpos": "start", "endpos": "end"})
            .with_columns(cn=pl.col("nMajor") + pl.col("nMinor"))
            .with_columns(pl.lit("ascat").alias("caller"))
        )
        return df


@define
class FreecResult(CopyNumberResult):
    schema: pa.DataFrameSchema = pa.DataFrameSchema(
        {
            "chr": pa.Column(pl.String),
            "start": pa.Column(pl.Int64),
            "end": pa.Column(pl.Int64),
            "copy number": pa.Column(pl.Int64),
            "WilcoxonRankSumTestPvalue": pa.Column(pl.Float64),
            "KolmogorovSmirnovPvalue": pa.Column(pl.Float64),
        },
        coerce=True,
    )

    def _read_segmental(self, file: Path) -> pl.LazyFrame:
        "Reads .p.value.txt"
        thresh: float = self.cfg.get("p_threshold", 0.05)
        df = (
            pl.scan_csv(file, separator="\t", infer_schema_length=None)
            .pipe(self.schema.validate)
            .filter(
                (pl.col("WilcoxonRankSumTestPvalue") <= thresh)
                & (pl.col("KolmogorovSmirnovPvalue") <= thresh)
            )
            .select(["chr", "start", "end", "copy number"])
            .rename({"copy number": "cn"})
            .with_columns(pl.lit("control_freec").alias("caller"))
        )
        return df


@define
class CNVKitResult(CopyNumberResult):
    schema: pa.DataFrameSchema = pa.DataFrameSchema(
        {
            "chromosome": pa.Column(pl.String),
            "start": pa.Column(pl.Int64),
            "end": pa.Column(pl.Int64),
            "cn": pa.Column(pl.Int64),
            "p_ttest": pa.Column(pl.Float64),
        },
        coerce=True,
    )

    def _read_segmental(self, file: Path) -> pl.LazyFrame:
        "Read and filter `call.cns` file"
        mappings: dict = Ensembl116.name_mapping(
            "ensembl_transcript", "gene_name", versioned=False
        )
        mappings.update(Ensembl116.name_mapping("mane_select", "gene_name"))
        lf: pl.LazyFrame = (
            pl.scan_csv(file, separator="\t", infer_schema_length=None)
            .pipe(self.schema.validate)
            .filter(pl.col("p_ttest") <= self.cfg.get("p_threshold", 0.05))
            .select(["chromosome", "start", "end", "gene", "cn"])
            .with_columns(pl.col("gene").str.split(","))
            .explode("gene")
            .filter(
                ~((pl.col("gene").str.starts_with("ccds|")) | (pl.col("gene") == "-"))
            )
            .with_columns(
                pl.when(pl.col("gene").str.starts_with("miRNA|"))
                .then(pl.col("gene").str.replace("miRNA|", ""))
                .otherwise(pl.col("gene").str.replace("^.*\\|", "").replace(mappings))
            )
            .group_by(["chromosome", "start", "end", "cn"])
            .agg(pl.col("gene"))
            .with_columns(
                pl.col("gene")
                .list.filter(~pl.element().str.contains("\\|"))
                .list.unique(),
                pl.lit("cnvkit").alias("caller"),
            )
            .rename({"chromosome": "chr", "gene": "genes"})
        )
        if (lf.collect()["genes"].list.len() == 0).all():
            lf = lf.drop("genes")
        else:
            lf = lf.with_columns(pl.col("genes").list.len().alias("n_genes"))
        return lf


def msisensor_pro_fmt(msisensor_path: Path) -> tuple:
    wanted_cols = ["Locus"] + list(Rename.repeat.values())
    df = pl.read_csv(
        msisensor_path, separator="\t", null_values="NA", infer_schema_length=None
    )
    if df.is_empty():
        return tuple([empty_table(Rename.repeat.values())] * 3)
    if "repeat_unit_bases" not in df.columns:
        df = df.with_columns(repeat_unit_bases=pl.lit(None))
    df = (
        add_loc(df, start_col="Start", end_col="End")
        .with_columns(
            pl.col("ClinGen_report").map_elements(
                get_clingen_link, return_dtype=pl.String
            )
        )
        .filter(pl.col("gene_name").is_not_null())
        .with_columns(
            pl.struct(s="source", acc="accession")
            .map_elements(
                lambda x: (
                    dbvar_link(x["acc"]) if x["acc"] != "NA" and x["acc"] else x["s"]
                ),
                return_dtype=pl.String,
            )
            .alias("source")
        )
        .rename(Rename.repeat)
        .with_columns(cs.by_dtype(pl.String).replace("NA", None))
        .cast(pl.String)
        .fill_null("-")
    )
    relevant = df.filter(pl.col("ClinGen") != "-").select(wanted_cols)
    nonrelevant = df.filter(pl.col("ClinGen") == "-").select(wanted_cols)
    return df, relevant, nonrelevant


def classify_cnv_fmt(
    path: Path,
    cnv_data_merge: Path,
    cols_keep: tuple[str] | None = ("cn", "callers"),
    chr_col: str = "chr",
    start_col: str = "start",
    end_col: str = "end",
) -> tuple:
    """
    Format cnv data from ClassifyCNV for the report

    Parameters
    ----------
    cnv_data_merge : Path | None
        Additional tabular data to add to ClassifyCNV output
    cols_keep : list[str] | None
        Columns from `cnv_data_merge` to retain
    """
    cnv = (
        pl.read_csv(path, separator="\t", null_values="NA", infer_schema_length=None)
        .filter(pl.col("Start").is_not_null() & pl.col("End").is_not_null())
        .with_columns(
            pl.struct(s="source", acc="accession")
            .map_elements(
                lambda x: (
                    dbvar_link(x["acc"]) if x["acc"] != "NA" and x["acc"] else x["s"]
                ),
                return_dtype=pl.String,
            )
            .alias("source"),
            pl.col("ClinGen_report").map_elements(
                get_clingen_link, return_dtype=pl.String
            ),
        )
    )

    def limit_genes(col, n):
        return pl.col(col).str.split(",").list.head(n).list.join(",")

    if cnv.is_empty():
        return tuple([empty_table(Rename.cnv)] * 3)

    # "group_by" operation required because data was exploded to become longer on
    #   dosage-sensitive genes in the R script
    wanted_cols = ["Locus"] + list(Rename.cnv.values()) + list(cols_keep)
    cnv = (
        add_loc(
            cnv,
            start_col="Start",
            end_col="End",
            chr_col="Chromosome",
        )
        .rename(Rename.cnv)
        .group_by("Locus")
        .agg((~cs.by_name("ClinGen")).unique().first(), pl.col("ClinGen"))
        .with_columns(
            pl.col("ClinGen").list.join(", "),
            limit_genes(Rename.cnv["All protein coding genes"], 25),
        )
    )
    to_merge = (
        pl.read_csv(cnv_data_merge)
        .pipe(add_loc, chr_col=chr_col, start_col=start_col, end_col=end_col)
        .select("Locus", *cols_keep)
    )
    with_cn = (
        cnv.join(to_merge, on="Locus", how="left")
        .select(wanted_cols)
        .cast(pl.String)
        .with_columns(cs.by_dtype(pl.String).str.replace("^$", "-"))
        .fill_null("-")
    )
    relevant = with_cn.filter(
        (pl.col("ClinGen") != "-") & (pl.col("ClinGen").str.len_chars() > 1)
    )
    nonrelevant = with_cn.filter(~pl.col("Locus").is_in(relevant["Locus"]))
    return with_cn, relevant, nonrelevant


def format_var_type(var: str) -> str:
    to_remove = ["_variant$"]
    for pat in to_remove:
        var = re.sub(pat, "", var)
    return var.strip()


def get_vcf_nth(col: str, n: int, cast_to: pl.DataType = pl.Int64) -> pl.Expr:
    if cast_to in {pl.Int64, pl.Float64}:
        return (
            pl.col(col)
            .str.split(",")
            .list.get(n, null_on_oob=True)
            .replace(".", None)
            .cast(cast_to)
        )
    return pl.col(str).str.split(",").list.get(n, null_on_oob=True).cast(cast_to)


def collapse_hgvs(
    df: pl.DataFrame | pl.LazyFrame,
    hgvsp: str = "HGVSp",
    hgvsc: str = "HGVSc",
    hgvsg: str = "HGVSg",
    col_added: str = "HGVS",
) -> pl.DataFrame:
    """
    Group HGVS columns (HGVSp, HGVSg, HGVSc) into a single column "HGVS"
    When available, protein-level and transcript-level identifiers are preferred over
    genomic
    """
    if isinstance(df, pl.DataFrame):
        ordering = [o for o in [hgvsp, hgvsc, hgvsg] if o in df.columns]
    else:
        ordering = [
            o for o in [hgvsp, hgvsc, hgvsg] if o in df.collect_schema().names()
        ]
    if not ordering:
        print("WARNING: no HGVS columns are in `df`")
        return df.with_columns(pl.lit(None).alias(col_added).cast(pl.String))
    return df.with_columns(pl.coalesce(ordering).alias(col_added).cast(pl.String))


def vep_fmt(
    vep_path: Path,
    variant_class: Literal["sv", "small"],
    format_fields: Path,
    sample: str | None = None,
    existing_variants: list[str] | None = None,
    mavedb: MaveDB | None = None,
    significant_hgvs: set | None = None,
) -> tuple:
    """Format and filter vep output into a dataframe with values ready to write
        into a reportlab table
    :param: existing_variants Accept only variants whose ids are in this list
        i.e. dbSNP, COSMIC
    """
    significant_hgvs = significant_hgvs or set()
    lf: pl.LazyFrame = pl.scan_csv(
        vep_path, separator="\t", null_values=["NA", "."], infer_schema_length=None
    )
    join_cols = ["CHROM", "POS", "ALT", "REF"]
    f_lf = pl.LazyFrame = pl.scan_csv(
        format_fields,
        separator="\t",
        null_values=["NA", "."],
        infer_schema_length=None,
    ).join(lf.select(join_cols), on=join_cols, how="right")
    format_cols = f_lf.collect_schema().names()
    if variant_class == "sv":
        rename: dict = deepcopy(Rename.sv)
    else:
        rename: dict = deepcopy(Rename.snp)
    if "sample" in format_cols and sample:
        f_lf = f_lf.filter(pl.col("sample") == sample)
    ad_cols = [col for col in format_cols if col.endswith("_AD")]
    af_cols = [col for col in format_cols if col.endswith("_AF")]
    if ad_cols:
        f_lf = f_lf.with_columns(
            get_vcf_nth(col, 1).alias(col) for col in ad_cols
        ).with_columns(
            pl.mean_horizontal(ad_cols).alias("Alt_depth"),
        )
    else:
        f_lf = f_lf.with_columns(pl.lit(None).alias("Alt_depth"))
    if af_cols:
        f_lf = f_lf.with_columns(pl.mean_horizontal(af_cols).alias("VAF"))
    else:
        f_lf = f_lf.with_columns(pl.lit(None).alias("VAF"))

    lf = lf.join(
        f_lf.select(join_cols + ["Alt_depth", "VAF"]), how="left", on=join_cols
    ).with_columns(
        pl.concat_str(["CHROM", "POS"], separator=":").alias("Loc"),
    )
    if lf.collect().is_empty():
        return tuple([empty_table(rename.values())] * 3)
    if existing_variants:
        lf = lf.filter(
            pl.col("Existing_variation")
            .str.split(";")
            .list.set_intersection(existing_variants)
            .len()
            >= 1
        )
    clinsig: str = Rename.snp["CLIN_SIG"]
    wanted_cols: list = list(rename.values())
    lf = lf.drop("Gene").pipe(collapse_hgvs).rename(rename)
    if variant_class != "sv" and ad_cols:
        lf = lf.with_columns(pl.col(rename["Alt_depth"]).cast(pl.Float64).round(0))
    if variant_class != "sv" and af_cols:
        lf = lf.with_columns((pl.col(rename["VAF"]).cast(pl.Float64) * 100).round(2))
    lf = (
        lf.with_columns(pl.col("HGVS").str.extract(r".*:(.*)$", 1))
        .pipe(lambda x: str_split_unique(x, "Variant Type", ";", ";", format_var_type))
        .with_columns(
            cs.by_dtype(pl.String)
            .str.replace_all("&", ", ", literal=True)
            .str.replace_all("_", " ", literal=True),
        )
        .with_columns(
            impact_score=pl.col("IMPACT").replace_strict(IMPACT_MAP, default=0)
        )
        .sort(pl.col("impact_score"), descending=True)
        .group_by(pl.col(["Gene", "Locus"]))
        .agg(pl.all().first())
        .filter((pl.col("Gene").is_not_null()) & (pl.col(clinsig) != "-"))
        .cast(pl.String)
        .fill_null("-")
    )
    if variant_class != "sv":
        lf = lf.filter(
            ~str_split_contains(clinsig, ["benign", "likely benign"], separator=";")
        )
    if not lf.collect().is_empty():
        df: pl.DataFrame = pl.concat(
            [
                sort_clinsig(d, clinsig_col=clinsig, separator=";")
                for d in lf.collect().partition_by("impact_score", maintain_order=True)
            ],
            how="vertical_relaxed",
        ).with_columns(cs.by_dtype(pl.String).str.replace_all(";", ", "))
        relevant_mask = df[clinsig].str.contains("pathogenic") | df["HGVS"].is_in(
            significant_hgvs
        )
        relevant = df.filter(relevant_mask)
        nonrelevant = df.filter(~relevant_mask)
        # if variant_class != "sv":
        #     nonrelevant = filter_multiallelic(
        #         nonrelevant,
        #         rename["VAF"],
        #         lambda x: all(list(map(lambda vaf: float(vaf) >= 0.6, x))),
        #         separator=",",
        #     )
        # if mavedb:
        #     nonrelevant = fmt_mavedb(nonrelevant, mavedb)
        return lf, relevant.select(wanted_cols), nonrelevant.select(wanted_cols)
    return tuple([empty_table(rename.values())] * 3)


def fmt_mavedb(
    df: pl.DataFrame,
    mavedb: MaveDB,
    join_spec: dict[str, str],
    gene_col: str = "Gene",
    mave_gene_col: str = "targetGenes.HGNC",
) -> pl.DataFrame:
    """
    Annotate `df` with columns from MaveDB lookups, ready to
    display in a reportlab table

    Parameters
    ----------
    hgvs_spec : dict
        Dictionary mapping columns in `df` to those in `mavedb` for joining.
        E.g. {"HGVS": "hgvs.g", "HGVS": "hgvs.c"}

    TODO: [2026-11-03 Tue] Finalize all the external variation sources
    before doing this. Also consider producing an HTML report
    for the uncertain variants or something similar.
    Would be way easier to manage

    Returns
    -------
    "[MaveDB ('SCORE_SET_TITLE')](SCORE_SET_URN): [FUNCTIONAL_SCORE](VARIANT_URN), "
    DataFrame with a
        - "Score"
        - "MaveDB Score Set"
        - "MaveDB classification"
    """
    tmp = []
    for gene in df[gene_col].unique():
        cur = df.filter(pl.col(gene_col) == gene)
        from_mavedb = mavedb.gene_score_set(gene, col=mave_gene_col)
        if from_mavedb.is_empty():
            tmp.append(cur)
            continue
        for left, right in join_spec.items():
            if right not in from_mavedb.columns:
                continue
            cleaned = (
                from_mavedb.filter(pl.col(right).is_not_null())
                .unique(right)
                .with_columns(pl.col(right).str.replace(".*:", ""))
            )
            cur = cur.join(
                cleaned,
                left_on=left,
                right_on=right,
                how="left",
            )
        tmp.append(cur)
    return pl.concat(tmp, how="diagonal_relaxed")


def str_split_contains(colname: str, query: list[str], separator: str = ";") -> pl.Expr:
    """Split the column `colname` by separator, and return a boolean series
    that is true if the row when split overlaps with query
    """
    return (
        pl.col(colname).str.split(separator).list.set_intersection(query).list.len()
    ) > 0


def sort_clinsig(
    df: pl.DataFrame, clinsig_col: str = "CLIN_SIG", separator: str = ""
) -> pl.DataFrame:
    if not separator:
        df = df.with_columns(
            clinsig_score=pl.col(clinsig_col).replace_strict(CLINSIG_MAP, default=0)
        )
    else:

        def find_max(terms: pl.Series) -> int:
            if terms.is_empty():
                return 0
            return max((map(lambda x: CLINSIG_MAP.get(x, 0), terms)))

        df = df.with_columns(
            clinsig_score=pl.col(clinsig_col)
            .str.split(separator)
            .map_elements(find_max, return_dtype=pl.Int64)
        )
    return df.sort(pl.col(clinsig_col), descending=True).drop("clinsig_score")


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


def empty_table(columns: Iterable) -> pl.DataFrame:
    schema = {c: pl.String for c in columns}
    return pl.DataFrame(schema=schema)


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
        current: pl.DataFrame = v["df"]
        if current.is_empty():
            continue
        current = current.select(gene_col)
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

    wanted_cols: tuple = TherapyDB.shared_cols + ("type",)
    final_cols = (
        "Therapy",
        "PubChemId",
        "Evidence category",
        "Evidence in sample",
        "Relevant cancers",
        "Study source",
        "Database source",
    )
    if len(dfs) > 0 and therapy_df.shape[0] > 0:
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
            .pipe(empty_string2null)
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
            .filter(pl.col("therapies").str.contains("[A-Za-z]+"))
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
                pl.col("therapies")
                .str.to_lowercase()
                .map_elements(str.capitalize, return_dtype=pl.String),
            )
            .rename(Rename.therapy)
            .select(list((Rename.therapy).values()))
            .cast(pl.String)
            .pipe(empty_string2null)
            .fill_null("-")
        )
        therapy_df = therapy_df.with_columns(
            cs.by_dtype(pl.String).str.replace_all(";", ", ")
        )
    else:
        therapy_df = pl.DataFrame({c: [] for c in final_cols})
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
        .select(wanted_cols)
        .rename({"Count": f"Count (Total: {total})"})
        .cast(pl.String)
        .fill_null("-")
    )
    return df


def sigprofiler_fmt(
    soln_activities_path: str,
    cosmic_signatures: str,
    abs_threshold: int = 100,
    rel_threshold: float = 0.25,
    excluded_signatures: str = "",
) -> list[pl.DataFrame]:
    """Format SigProfilerAssignment results for reportlab

    :param soln_activities_path: Path to SigProfiler "Assignment_Solution_Activities.txt" file
    :param cosmic_signatures: Path to csv file with COSMIC signature metadata

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
    df: pl.DataFrame = pl.read_csv(
        soln_activities_path, separator="\t", infer_schema_length=None
    ).select(~cs.by_name(*excluded))
    signatures: list = df.columns[1:]
    if df.is_empty():
        return []
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
        df.with_columns((sig_cols > abs_threshold) | ((sig_cols / m) > rel_threshold))
        .with_columns(replace_expr)
        .with_columns(Signatures=pl.concat_list(signatures).list.drop_nulls())
        .select(["Samples", "Signatures"])
        .join(sums, on="Samples")
        .explode("Signatures")
        .join(frequencies, on=["Samples", "Signatures"])
        .join(counts, on=["Samples", "Signatures"])
        .group_by("Samples")
        .agg(pl.col(["Signatures", "Frequency", "Count"]), pl.col("m").first())
    )
    cosmic_data: pl.DataFrame = pl.read_csv(cosmic_signatures)
    formatted = [make_signature_table(s, cosmic_data) for s in filtered.iter_slices(1)]
    return formatted
