#!/usr/bin/env ipython

import os
import re
from io import StringIO
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

import polars as pl
from snakemake.script import snakemake as smk


def bcftools_stats(input, output) -> None:
    command = f"bcftools norm -a --atom-overlaps '*' -d none {input}"
    run(f"{command} | bcftools stats > {output}", shell=True)


def get_substitution_types(
    bcftools_stats_file: str, by_convention: bool = True
) -> pl.DataFrame:
    with open(bcftools_stats_file, "r") as f:
        text = f.read().splitlines()
    nucl_map: dict[str, str] = {"G": "C", "A": "T", "C": "G", "T": "A"}
    wanted_lines = []
    for index, line in enumerate(text):
        if (
            line == "# ST, Substitution types:"
            and text[index + 1] == "# ST	[2]id	[3]type	[4]count"
        ):
            i: int = index + 2
            while not text[i].startswith("#"):
                wanted_lines.append(i)
                i += 1
            break
    data = "\n".join(text[wanted_lines[0] : wanted_lines[-1] + 1])
    df = pl.read_csv(
        StringIO(data),
        has_header=False,
        separator="\t",
        new_columns=["ST", "ID", "type", "count"],
    )
    if by_convention:

        def convert(sub: str) -> str:
            if sub.startswith("C") or sub.startswith("T"):
                return sub
            return ">".join(list(map(lambda x: nucl_map.get(x), sub.split(">"))))

        df = (
            df.with_columns(
                pl.col("type").map_elements(convert, return_dtype=pl.String)
            )
            .group_by("type")
            .agg(pl.col("count").sum(), pl.col("ST", "ID").first())
            .select(df.columns)
            .sort("type")
        )
    return df


data: Path = Path(smk.params["data_path"])
target = str(smk.output)
cases = smk.params["cases"]

dfs: list[pl.DataFrame] = []
cwd = os.getcwd()

found_cases = []
files = []
for case in cases:
    case_dir = data.joinpath(case).resolve().absolute()
    look_for = next(case_dir.rglob("7-*-Small_high_conf.vcf.gz"), None)
    if look_for:
        files.append(look_for)
        found_cases.append(case)

with TemporaryDirectory() as tmpdir:
    os.chdir(tmpdir)
    for f, c in zip(files, cases):
        # sample_name = re.findall(r"7-(.*)-Small_high_conf.vcf", f.stem)[0]
        bcftools_stats(f, "tmp_stats.txt")
        df = get_substitution_types("tmp_stats.txt").with_columns(sample=pl.lit(c))
        dfs.append(df)
os.chdir(cwd)
final = pl.concat(dfs)
final.write_csv(target, separator="\t")
