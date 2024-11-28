#!/usr/bin/env ipython


from tempfile import TemporaryFile

import click
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

# rx method is equivalent to [], rx2 is [[]]


def r2pd(robject):
    with (ro.default_converter + pandas2ri.converter).context():
        return ro.conversion.get_conversion().rpy2py(robject)


@click.command()
@click.option("-i", "--input", required=True, help="Input file")
@click.option("-o", "--output", required=True, help="Output")
@click.option("-c", "--caller", required=True, help="Tool producing CNV calls")
def classify_cnv_format(caller: str, input: str, output: str):
    _classify_cnv(caller, input, output)


def _classify_cnv(caller: str, input: str, output: str, tumor_sample: str = ""):
    import polars as pl

    def dup_or_del(x):
        return "DEL" if x < 2 else "DUP"

    if caller.lower() == "cnvkit":
        df = pl.read_csv(input, separator="\t").select(
            ["chromosome", "start", "end", "cn"]
        )
    elif caller.lower() == "facets":
        base = importr("base")
        df: pl.DataFrame = (
            pl.from_pandas(r2pd(base.readRDS(input).rx2("segs")))
            .select(["chrom", "start", "end", "tcn.em"])
            .rename({"tcn.em": "cn"})
        )
    elif caller.lower() == "dellycnv":
        from subprocess import run

        if not tumor_sample:
            raise ValueError("Tumor sample name must be provided")
        qstring = f"bcftools query -f '%CHROM\t%POS0\t%END\t[%CN]\n' -s {tumor_sample} {input} > tmp.bed"
        print(qstring)
        run(
            qstring,
            shell=True,
        )
        df = pl.read_csv(
            "tmp.bed",
            separator="\t",
            new_columns=["chromosome", "start", "end", "cn"],
            dtypes=[pl.String, pl.Int64, pl.Int64, pl.Int64],
        )
        run("rm tmp.bed", shell=True)
    else:
        raise ValueError("Unknown caller provided")

    df.filter(pl.col("cn") != 2).with_columns(
        type=pl.col("cn").map_elements(dup_or_del, return_dtype=pl.String)
    ).drop("cn").write_csv(output, separator="\t", include_header=False)


@click.command()
@click.option("-n", "--normal", required=True)
@click.option("-t", "--tumor", required=True)
@click.option("-s", "--snps", help="Path to snps file", required=True)
@click.option(
    "-r", "--chr_files", help="Path to directory with chromosomes", required=True
)
@click.option("-l", "--chr_len_file", required=True)
@click.option("-o", "--output", required=True)
@click.option("-b", "--breakpoint_threshold", default=1.2)
@click.option("-m", "--minimal_coverage_per_position", default=0)
@click.option("-p", "--ploidy", default="2,3")
@click.option(
    "-r",
    "--read_type",
    default="paired_end",
    type=click.Choice(["paired_end", "mate_pair", "single"]),
)
@click.option("-w", "--window", default=0)
@click.option("-d", "--outputdir", default=".")
@click.option("-c", "--contamination", default=0)
@click.option("-e", "--sex", default="XY")
@click.option("-y", "--breakpoint_type", default=2)
@click.option("-x", "--max_threads", default=1)
@click.option("--noisy_data", default=True, type=bool, is_flag=True)
@click.option(
    "-f",
    "--file_format",
    default="BAM",
    help="input format of tumor and normal files",
)
@click.option("--print_na", default=False, type=bool, is_flag=True)
@click.option("-i", "--intervals", help="Exome target intervals")
def make_freec_config(
    tumor: str,
    normal: str,
    chr_len_file: str,
    chr_files: str,
    max_threads: int,
    minimal_coverage_per_position,
    window,
    file_format,
    ploidy,
    outputdir,
    contamination,
    breakpoint_threshold,
    breakpoint_type,
    sex,
    snps,
    noisy_data,
    print_na,
    read_type,
    intervals,
    output,
):
    import toml

    template: dict = {
        "general": {
            "chrLenFile": chr_len_file,
            "window": window,
            "ploidy": ploidy,
            "outputDir": outputdir,
            "sex": sex,
            "breakPointType": breakpoint_type,
            "contamination": contamination,
            "chrFiles": chr_files,
            "maxThreads": max_threads,
            "breakPointThreshold": breakpoint_threshold,
        },
        "sample": {"mateFile": normal, "inputFormat": file_format},
        "control": {"mateFile": tumor, "inputFormat": file_format},
        "BAF": {
            "SNPfile": snps,
            "minimalCoveragePerPosition": minimal_coverage_per_position,
        },
    }
    if read_type == "paired_end":
        template["sample"]["mateOrientation"] = "FR"
        template["control"]["mateOrientation"] = "FR"
    elif read_type == "mate_pair":
        template["sample"]["mateOrientation"] = "RF"
        template["control"]["mateOrientation"] = "RF"
    elif read_type == "single":
        template["sample"]["mateOrientation"] = 0
        template["control"]["mateOrientation"] = 0
    if noisy_data:
        template["general"]["noisyData"] = "TRUE"
    if not print_na:
        template["general"]["printNA"] = "false"
    if intervals:
        template["target"] = {"captureRegions": intervals}
        template["general"]["window"] = 0
    with TemporaryFile("r+") as w:
        toml.dump(template, w)
        w.seek(0)
        lines = w.readlines()
    with open(output, "w") as w:
        modified: list[str] = [l.replace('"', "").replace("'", "") for l in lines]
        w.write("".join(modified))
