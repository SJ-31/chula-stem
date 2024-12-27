import click
import polars as pl

from chula_stem.utils import empty_string2null, contain_join

IMPACT_MAP: dict = {
    "HIGH": 3,
    "MODERATE": 2,
    "LOW": 1,
    "MODIFIER": 1,
    None: 0,
}


def merge_variant_calls(
    df: pl.DataFrame,
    grouping_cols: list,
    tool_source_tag: str = "TOOL_SOURCE",
    minimum_callers: int = 2,
    vaf_adaptive: bool = False,
    separator: str = ";",
) -> pl.DataFrame:
    """Merge variant calling results

    By default, merge variant calling results with a "majority vote" strategy, where
    variants are accepted if they have been called by n = `minimum_callers` callers

    With `vaf_adaptive` mode, proposed by Wang et. al 2020 in SomaticCombiner
    variants are accepted if
    - called by Strelka + Mutect2 and 0.03 <= tumor VAF <= 0.1
    - called by Mutect2 and VAF < 0.03
    - called by >= minimum_callers

    :param: minimum_callers Variants must be called by this number of different callers
    :param: tool_source_tag col containing the variant caller names
    :param: canonical_only When multiple transcripts are available, choose only the
        canonical variant if it is present
    :param: highest_impact Choose the variants with the highest impact class

    :returns: filtered tsv file
    """
    original_cols: list = df.columns
    to_average: list = []
    if "VAF" in original_cols:
        to_average.append("VAF")
    if "Alt_depth" in original_cols:
        to_average.append("Alt_depth")
    vep_cols = list(
        filter(lambda x: not (x in grouping_cols or x in to_average), original_cols)
    )
    split_unique: list = vep_cols
    keep_all: list = to_average + split_unique
    unique_expr: list = [pl.col(u).list.unique() for u in split_unique]
    avg_expr: list = [pl.col(x).list.mean() for x in to_average]
    grouped = (
        df.group_by(grouping_cols)
        .agg(keep_all)
        .with_columns(avg_expr + unique_expr)
        .with_columns(n_callers=pl.col(tool_source_tag).list.len())
    )
    if vaf_adaptive:
        mutect_strelka_cols = [
            (pl.col(tool_source_tag).list.contains(c)).alias(f"has_{c}")
            for c in ["mutect2", "strelka2"]
        ]
        grouped = grouped.with_columns(mutect_strelka_cols).filter(
            (pl.col("has_mutect2") & (pl.col("VAF") < 0.03))
            | (
                ((pl.col("has_mutect2") & pl.col("has_strelka2")))
                & (pl.col("VAF") <= 0.1)
                & (pl.col("VAF") >= 0.03)
            )
            | (pl.col("n_callers") >= minimum_callers)
        )
    else:
        grouped = grouped.filter(pl.col("n_callers") >= minimum_callers)
    grouped = grouped.with_columns(
        pl.col(split_unique).cast(pl.List(pl.String)).list.join(separator)
    ).select(original_cols)
    return grouped


def resolve_transcripts(
    df: pl.DataFrame,
    grouping_cols: list,
    impact: bool = True,
    canonical: bool = True,
    informative: bool = True,
) -> pl.DataFrame:
    """Choose between transcripts of gene based on...
    - Highest impact
    - Whether or not it is the canonical transcript
    - Which is most informative/Whichever has the most non-empty cells
    """

    def by_impact(df: pl.DataFrame) -> pl.DataFrame:
        v = "IMPACT_VAL"
        df = (
            df.with_columns(pl.col("IMPACT").replace_strict(IMPACT_MAP).alias(v))
            .filter(pl.col(v) == pl.col(v).max())
            .drop(v)
        )
        return df

    def by_canonical(df: pl.DataFrame) -> pl.DataFrame:
        filtered = df.filter(pl.col("CANONICAL") == "YES")
        if filtered.is_empty():
            return df
        return filtered

    def by_informative(df: pl.DataFrame) -> pl.DataFrame:
        original_cols: list = df.columns
        sum_col: str = "na_count"
        df = (
            df.with_columns(pl.sum_horizontal(pl.all().is_null()).alias(sum_col))
            .filter(pl.col(sum_col) == pl.col(sum_col).min())
            .select(original_cols)
        )
        return df

    dfs: list[pl.DataFrame] = []
    for _, group in df.group_by(grouping_cols):
        group: pl.DataFrame
        if canonical and "CANONICAL" in df.columns:
            group = by_canonical(group)
        if impact:
            group = by_impact(group)
        if informative:
            group = by_informative(group)
        dfs.append(group)
    if dfs:
        resolved = pl.concat(dfs)
        return resolved
    return pl.DataFrame(schema=df.schema)


def standard_filters(
    df: pl.DataFrame,
    min_tumor_depth: int,
    max_normal_depth: int,
    min_vaf: float,
    accepted_filters: str,
) -> pl.DataFrame:
    filters = accepted_filters.split(",")
    if "Alt_depth_normal" in df.columns:
        df = df.cast({"Alt_depth_normal": pl.Float64}).filter(
            (pl.col("Alt_depth_normal") <= max_normal_depth)
            | (pl.col("Alt_depth_normal").is_null())
        )
    if "Alt_depth" in df.columns:
        df = df.cast({"Alt_depth": pl.Float64}).filter(
            (pl.col("Alt_depth") >= min_tumor_depth) | (pl.col("Alt_depth").is_null())
        )
    if "VAF" in df.columns:
        df = df.cast({"VAF": pl.Float64}).filter(
            (pl.col("VAF") >= min_vaf) | (pl.col("VAF").is_null())
        )
    if filters:
        df = df.filter(
            pl.col("FILTER").str.split(";").list.set_intersection(filters).list.len()
            >= 1
        )
    return df


def region_filter(
    df: pl.DataFrame,
    ignore_file: str,
    chr_col: int = 0,
    start_col: int = 1,
    end_col: int = 2,
) -> pl.DataFrame:
    """Remove variants within the specified regions

    Any variants on the same chromosome of `ignore_file` and within the range specified
        by `start_col` and `end_col` are filtered out

    :param ignore_file: a TSV file (e.g. assumed to be a BED file by default)
            containing regions to ignore
    """
    wanted_cols: list = ["chr", "start", "end"]
    original_cols = df.columns
    df = (
        df.with_columns(
            pl.col("Loc")
            .str.split(":")
            .list.to_struct(fields=["chr", "start"])
            .alias("chr_start")
        )
        .unnest("chr_start")
        .cast({"chr": pl.String, "start": pl.Int64})
    )
    ignore = (
        pl.read_csv(
            ignore_file, separator="\t", null_values="NA", infer_schema_length=None
        )
        .with_columns(
            pl.nth(chr_col).alias(wanted_cols[0]),
            pl.nth(start_col).alias(wanted_cols[1]),
            pl.nth(end_col).alias(wanted_cols[2]),
        )
        .select(wanted_cols)
        .cast({"chr": pl.String, "start": pl.Int64, "end": pl.Int64})
    )
    joined = contain_join(
        df, ignore, x_start="start", y_start="start_right", y_end="end", on=["chr"]
    )
    df = df.filter(~pl.col("Loc").is_in(joined["Loc"]))
    return df.select(original_cols)


@click.command()
@click.option("-o", "--output", required=True, help="Output tsv path")
@click.option("-i", "--input_tsv", required=True, help="Input tsv path")
@click.option("--min_tumor_depth", required=False, default=0)
@click.option("--max_normal_depth", required=False, default=50)
@click.option("--min_vaf", required=False, default=0)
@click.option("--impact", required=False, default=False, is_flag=True)
@click.option("--canonical", required=False, default=False, is_flag=True)
@click.option("--accepted_filters", required=False, default="PASS")
@click.option("--informative", required=False, default=False, is_flag=True)
@click.option("--min_callers", required=False, default=2)
@click.option("--vaf_adaptive", required=False, default=False, is_flag=True)
@click.option("--tool_source_tag", required=False, default="TOOL_SOURCE")
@click.option("-g", "--ignore_regions", required=False, default="")
@click.option(
    "--chr_col",
    required=False,
    default=0,
    help="Chromosome column index for `ignore_regions` file",
)
@click.option(
    "--start_col",
    required=False,
    default=1,
    help="Starting column index for `ignore_regions` file",
)
@click.option(
    "--end_col",
    required=False,
    default=2,
    help="Ending column index for `ignore_regions` file",
)
def qc_main(
    input_tsv: str,
    output: str,
    min_tumor_depth: int,
    max_normal_depth: int,
    min_vaf: float,
    accepted_filters: str,
    impact: bool,
    canonical: bool,
    informative: bool,
    min_callers: bool,
    vaf_adaptive: bool,
    tool_source_tag: str = "TOOL_SOURCE",
    ignore_regions: str = "",
    chr_col: int = 0,
    start_col: int = 1,
    end_col: int = 2,
):
    _qc_main(
        input_tsv=input_tsv,
        output=output,
        min_tumor_depth=min_tumor_depth,
        max_normal_depth=max_normal_depth,
        min_vaf=min_vaf,
        accepted_filters=accepted_filters,
        impact=impact,
        canonical=canonical,
        informative=informative,
        min_callers=min_callers,
        vaf_adaptive=vaf_adaptive,
        tool_source_tag=tool_source_tag,
        ignore_regions=ignore_regions,
        chr_col=chr_col,
        start_col=start_col,
        end_col=end_col,
    )


def _qc_main(
    input_tsv: str,
    output: str,
    min_tumor_depth: int,
    max_normal_depth: int,
    min_vaf: float,
    accepted_filters: str,
    impact: bool,
    canonical: bool,
    informative: bool,
    min_callers: bool,
    vaf_adaptive: bool,
    tool_source_tag: str = "TOOL_SOURCE",
    ignore_regions: str = "",
    chr_col: int = 0,
    start_col: int = 1,
    end_col: int = 2,
) -> None:
    df: pl.DataFrame = pl.read_csv(
        input_tsv,
        separator="\t",
        infer_schema_length=None,
        null_values=["NA", "."],
    )
    print(f"Original shape: {df.shape}")
    grouping_cols: list = ["Loc", "Ref", "Alt"]
    df = standard_filters(
        df,
        min_tumor_depth=min_tumor_depth,
        max_normal_depth=max_normal_depth,
        min_vaf=min_vaf,
        accepted_filters=accepted_filters,
    )
    print(f"After standard_filters: {df.shape}")
    df = resolve_transcripts(
        df,
        grouping_cols,
        impact=impact,
        canonical=canonical,
        informative=informative,
    )
    print(f"After resolve_transcripts: {df.shape}")
    df = merge_variant_calls(
        df,
        grouping_cols,
        tool_source_tag=tool_source_tag,
        minimum_callers=min_callers,
        vaf_adaptive=vaf_adaptive,
    )
    print(f"After merge_variant_calls: {df.shape}")
    if ignore_regions:
        print(f"Before filtering by regions in {ignore_regions}: {df.shape}")
        df = region_filter(df, ignore_regions, chr_col, start_col, end_col)
        print(f"After filtering by regions in {ignore_regions}: {df.shape}")
    df.unique(["Loc", "Feature"], keep="first", maintain_order=True).pipe(
        empty_string2null
    ).write_csv(output, separator="\t", null_value="NA")
