import os
import re
from io import StringIO
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

import polars as pl

outpath: Path = Path("/home/shannc/Bio_SDD/chula-stem/analyses/data_all/output/PDAC")
target = "/home/shannc/Bio_SDD/chula-stem/analyses/output/pdac_sbs.tsv"
target2 = "/home/shannc/Bio_SDD/chula-stem/analyses/output/pdac_sbs_all.tsv"


def bcftools_stats(input, output) -> None:
    run(f"bcftools stats {input} > {output}", shell=True)


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


dfs: list[pl.DataFrame] = []
dfs2: list[pl.DataFrame] = []
cwd = os.getcwd()
if outpath.exists():
    files = list(outpath.rglob("7-P*-Small_high_conf.vcf.gz"))
    with TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        for f in files:
            sample_name = re.findall(r"-(P[0-9_]+)-", f.stem)[0]
            bcftools_stats(f, "tmp_stats.txt")
            df = get_substitution_types("tmp_stats.txt").with_columns(
                sample=pl.lit(sample_name)
            )
            df2 = get_substitution_types("tmp_stats.txt", False).with_columns(
                sample=pl.lit(sample_name)
            )
            dfs.append(df)
            dfs2.append(df2)
    os.chdir(cwd)
    final = pl.concat(dfs)
    final.write_csv(target, separator="\t")
    final2 = pl.concat(dfs2)
    final2.write_csv(target2, separator="\t")

old_path: Path = Path("/data/project/stemcell/PDAC/Exome_PDAC_HN00215149")
if old_path.exists():
    old_dfs = []
    with TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        for f in old_path.rglob("*final.vcf"):
            sample_name = f.stem.replace(".final.vcf", "")
            bcftools_stats(f, "tmp_stats.txt")
            odf = get_substitution_types("tmp_stats.txt").with_columns(
                sample=pl.lit(sample_name)
            )
            old_dfs.append(odf)
    os.chdir(cwd)
    ptarget = "/data/home/shannc/chula-stem/analyses/output/pdac_sbs_previous.tsv"
    old_final = pl.concat(old_dfs)
    old_final.write_csv(ptarget)
