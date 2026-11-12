#!/usr/bin/env ipython

import argparse
from pathlib import Path

import great_tables as gt
import mudata as md
import polars as pl


def format_df(
    df: pl.DataFrame, n_top: int = 10, excluded: list | None = None
) -> pl.DataFrame:
    excluded = excluded or []
    get_ordered_exprs = [
        pl.col("clones").list.eval(pl.element().struct.field(e)).alias(e)
        for e in ["clone_id", "size", "receptor_type", "receptor_subtype"]
    ]
    return (
        df.filter(~pl.col("Sample_Name").is_in(excluded))
        .with_columns(
            pl.struct(
                clone_id=pl.col("clone_id"),
                size=pl.len().over("Sample_Name", "clone_id"),
                receptor_type=pl.col("receptor_type").cast(pl.String),
                receptor_subtype=pl.col("receptor_subtype").cast(pl.String),
            ).alias("clones")
        )
        .group_by("Sample_Name")
        .agg(pl.col("clone_id"), pl.col("clones"))
        .with_columns(
            pl.col("clones")
            .list.eval(pl.element().struct.field("size"))
            .list.sum()
            .alias("clone_total")
        )
        .with_columns(
            pl.col("clone_id").list.unique().list.len().alias("n_clones"),
            pl.col("clones")
            .list.unique()
            .list.eval(
                pl.element().sort_by(pl.element().struct.field("size"), descending=True)
            )
            .list.head(n_top),
        )
        .drop("clone_id")
        .with_columns(*get_ordered_exprs)
        .cast({"size": pl.Array(pl.Float64, n_top)})
        .with_columns(
            (pl.col("size") / pl.col("clone_total")).alias("proportion"),
            pl.col("size").alias("size_raw"),
        )
        .with_columns(
            pl.col("size_raw")
            .arr.eval(pl.element().cast(pl.Int64).cast(pl.String))
            .arr.join(", "),
            pl.when(pl.col("receptor_type").list.n_unique() == 1)
            .then("All " + pl.col("receptor_type").list.unique().list.first())
            .otherwise(pl.col("receptor_type").list.join(", "))
            .alias("receptor_type"),
            pl.when(pl.col("receptor_subtype").list.n_unique() == 1)
            .then("All " + pl.col("receptor_subtype").list.unique().list.first())
            .otherwise(pl.col("receptor_subtype").list.join(", "))
            .alias("receptor_subtype"),
        )
    )


def make_table(df: pl.DataFrame, n_top: int = 10) -> gt.GT:
    cols_center = ["size_raw", "receptor_type", "receptor_subtype"]
    table = (
        gt.GT(df, groupname_col="Sample_Name")
        .tab_header("Cohort Report")
        .cols_hide(["clones", "clone_id"])
        .tab_source_note(
            f"Values are for the top {n_top} largest clones in each sample, ordered"
        )
        .tab_footnote(
            "The dashed horizontal bar is the maximum proportion",
            locations=gt.loc.column_labels("proportion"),
        )
        .cols_align("center", cols_center)
        .fmt_integer(columns=["clone_total", "n_clones"])
        .fmt_nanoplot("proportion", plot_type="bar", reference_line=1, alt="Proportion")
        .fmt_nanoplot("size", plot_type="bar", autoscale=True)
        .cols_move_to_start(["size", "proportion"])
        .opt_stylize(style=6)
        .cols_label(
            size="Size comparison",
            size_raw="Sizes",
            clone_total="Total clones",
            proportion="Proportion of total",
            n_clones="N. unique clones",
            receptor_type="Receptor types",
            receptor_subtype="Receptor subtypes",
        )
    )
    return table


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument(
        "-n",
        "--n_top",
        default=10,
        help="Number of top clones per sample to report",
        action="store",
        type=int,
    )
    parser.add_argument("-o", "--output", required=True, type=Path)
    parser.add_argument(
        "-s",
        "--sample_col",
        default="Sample_Name",
        help="Sample name column in mudata object",
        action="store",
    )
    parser.add_argument(
        "-a",
        "--airr_mod",
        default="airr",
        help="Slot in mudata object containing AIRR data",
        action="store",
    )
    parser.add_argument(
        "-e",
        "--excluded",
        default="Undetermined,Multiplet",
        help="Comma-delimited list of samples to exclude",
        action="store",
    )
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__":
    args = parse_args()
    mdata = md.read_h5mu(args["input"])
    to_exclude = args["excluded"].split(",") if args["excluded"] else []
    df: pl.DataFrame = pl.from_pandas(mdata[args["airr_mod"]].obs.reset_index())
    formatted = format_df(df, n_top=args["n_top"], excluded=to_exclude)
    table = make_table(formatted, n_top=args["n_top"])
    output: Path = args["output"]
    if output.suffix == ".html":
        table.write_raw_html(output, inline_css=True)
    else:
        table.gtsave(output)
