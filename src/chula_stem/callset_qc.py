import polars as pl
import polars.selectors as cs


def print_df(df: pl.DataFrame) -> pl.DataFrame:
    print(df)
    return df


def save_df(df: pl.DataFrame) -> pl.DataFrame:
    """Save the current polars dataframe to a temporary file for debugging purposes
    Can be used in a pipe
    """
    import datetime
    import os
    import sys

    import polars as pl
    import polars.selectors as cs

    time = datetime.datetime.now().strftime("%Y-%m-%d-%M_%S")
    out = f"{time}_{os.path.basename(sys.argv[0])}"
    try:
        df.write_csv(f"{out}.tsv", separator="\t")
    except pl.exceptions.ComputeError:
        df.write_json(f"{out}.json", pretty=True)
        os.unlink(f"{out}.tsv")
    return df


def merge_variant_calls(
    df: pl.DataFrame,
    grouping_cols: list,
    tool_source_tag: str = "TOOL_SOURCE",
    minimum_callers: int = 3,
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
    original_shape: tuple = df.shape
    original_cols: list = df.columns
    to_average: list = ["VAF", "Alt_depth"]
    vep_cols = list(
        filter(
            lambda x: not (x in grouping_cols or x in to_average), original_cols
        )
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
            (pl.col(tool_source_tag) == c).alias(f"has_{c}")
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
        pl.col(split_unique).list.join(separator)
    ).select(original_cols)
    new_shape = grouped.shape
    print(f"Shape before merging: {original_shape}\nShape after: {new_shape}")
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
        impact_map: dict = {
            "HIGH": 3,
            "MODERATE": 2,
            "LOW": 1,
            "MODIFIER": 1,
            None: 0,
        }
        v = "IMPACT_VAL"
        df = (
            df.with_columns(pl.col("IMPACT").replace_strict(impact_map).alias(v))
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
        str_cols: list = df.select(cs.by_dtype(pl.String)).columns
        pref: str = "__notna__"
        replace_expr = [
            pl.col(s).replace_strict({"NA": 0}, default=1).alias(f"{pref}{s}")
            for s in str_cols
        ]
        sum_col: str = "notna_sum"
        df = (
            df.with_columns(replace_expr)
            .with_columns(pl.sum_horizontal(cs.starts_with(pref)).alias(sum_col))
            .filter(pl.col(sum_col) == pl.col(sum_col).max())
            .select(original_cols)
        )
        return df

    dfs: list[pl.DataFrame] = []
    for group in df.group_by(grouping_cols):
        group: pl.DataFrame
        if canonical and "CANONICAL" in df.columns:
            group = by_canonical(group)
        if impact:
            group = by_impact(group)
        if informative:
            group = by_informative(group)
        dfs.append(group)
    resolved = pl.concat(dfs)
    return resolved


def qc_main(input_tsv: str) -> None:
    df: pl.DataFrame = pl.read_csv(input_tsv, separator="\t")
    grouping_cols: list = ["Loc", "Ref", "Alt"]
    df = merge_variant_calls(df, grouping_cols)
    df = resolve_transcripts(df, grouping_cols)
