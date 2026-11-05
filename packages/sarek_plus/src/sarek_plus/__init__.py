#!/usr/bin/env python3

import copy
import json
from itertools import product
from pathlib import Path
from typing import Literal, TypeAlias

import click
import yaml
from attrs import Factory, define, field

from sarek_plus.databases import get_therapy_df
from sarek_plus.report.main import VariantCallingReport

OUTPUT_TYPE: TypeAlias = Literal["tumor_only", "paired", "germline"]


@click.command()
@click.option("-o", "--output", required=False, help="Name of output pdf")
@click.option(
    "-t",
    "--report_type",
    required=True,
    default=None,
    help="Type of report to be generated",
)
@click.option(
    "-s", "--specification", required=True, help="Path specification file (YAML)"
)
@click.option(
    "-d",
    "--tmpdir",
    required=False,
    default="./report_tmp",
    help="Path to temporary directory for report files",
)
def entry_point(report_type: str, specification: str, output: str, tmpdir: str) -> None:
    with open(specification, "r") as f:
        spec: dict = json.load(f)
    try:
        paths: dict = spec["paths"]
        assert isinstance(paths, dict)
    except (KeyError, AssertionError):
        raise ValueError("`paths` map must be specified!")
    out = output or paths["output"]
    kws = spec.get("misc", {})
    report = None
    if report_type == "variant_calling":
        report = VariantCallingReport(
            filename=out,
            metadata=spec.get("meta", {}),
            other_text=spec.get("other_text", {}),
            tmpdir=tmpdir,
            paths=spec["paths"],
            **kws,
        )
    if report:
        report.build()
        report.merge()


@click.command()
@click.option(
    "-i",
    "--input-files",
    required=True,
    help="Path to directory containing input files",
)
@click.option(
    "-c",
    "--gene-col",
    default="Gene",
    show_default=True,
    help="For dataframe input, the column containing gene names",
)
@click.option(
    "-f",
    "--filter-confident",
    is_flag=True,
    default=False,
    show_default=True,
    help="Whether or not to filter confident therapy entries",
)
@click.option(
    "-n",
    "--not-in-db",
    default="not_in_databases.txt",
    show_default=True,
    help="File to write genes that weren't found",
)
@click.option(
    "-i",
    "--ignore_failed_file",
    required=True,
    help="If supplied, will ignore all genes in this file",
)
@click.option(
    "-e", "--cache", required=True, help="DB file to cache data", default=None
)
@click.option(
    "-d",
    "--databases",
    required=True,
    help="Comma-delimited string specifying which databases to include in the lookup",
    default=None,
    type=str,
)
def therapy_cache_from_files(
    input_files: str,
    not_in_db: str,
    filter_confident: bool,
    cache: str,
    gene_col: str = "Genes",
):
    import polars as pl

    genes: set = set()
    for path in Path(input_files).rglob("*"):
        path: Path
        if path.is_dir():
            continue
        if path.suffix == ".tsv":
            df = pl.read_csv(path, separator="\t", infer_schema_length=None)
            genes |= set(df[gene_col])
        elif path.suffix == ".csv":
            df = pl.read_csv(path, infer_schema_length=None)
            genes |= set(df[gene_col])
        else:
            with open(path, "r") as f:
                genes |= set(f.read().splitlines())
    _ = get_therapy_df(
        genes,
        cache=Path(cache),
        not_in_db_file=not_in_db,
        filter_confident=filter_confident,
        lookup=True,
    )


@click.command()
@click.option("-i", "--input", required=True, help="Input file")
@click.option("-o", "--output", required=True, help="Output")
@click.option("-c", "--caller", required=True, help="Tool producing CNV calls")
@click.option(
    "-t", "--tumor_sample", required=False, help="Name of tumor sample", default=""
)
def classify_cnv_format(caller: str, input: str, output: str, tumor_sample: str):
    import polars as pl

    def dup_or_del(x):
        return "DEL" if x < 2 else "DUP"

    if caller.lower() == "cnvkit":
        df = pl.read_csv(input, separator="\t", infer_schema_length=None).select(
            ["chromosome", "start", "end", "cn"]
        )
    elif caller.lower() == "dellycnv":
        from subprocess import run

        if not tumor_sample:
            raise ValueError("Tumor sample name must be provided")
        qstring = f"bcftools query -f '%CHROM\t%POS0\t%END\t[%CN]\n' -s {tumor_sample} {input} > tmp.bed"
        run(
            qstring,
            shell=True,
        )
        df = pl.read_csv(
            "tmp.bed", separator="\t", new_columns=["chromosome", "start", "end", "cn"]
        )
        run("rm tmp.bed", shell=True)
    else:
        raise ValueError("Unknown caller provided")

    df.filter(pl.col("cn") != 2).with_columns(
        type=pl.col("cn").map_elements(dup_or_del, return_dtype=pl.String)
    ).drop("cn").write_csv(output, separator="\t", include_header=False)


def group_samples(samples: list[dict], config: dict) -> list[tuple[str, str]]:
    """
    Group `samples` by status and return a list of (normal, tumor) pairs.
    If multiple tumor and normal samples are present, they are crossed
    with `itertools.product`. If no samples of a given status are present,
    the empty string is left instead.
    """

    normals, tumors = [], []
    for s in samples:
        if s["status"] == "0":
            normals.append(s["sample"])
        else:
            tumors.append(s["sample"])
    if len(normals) == 0:
        return [("", t) for t in tumors]
    elif len(tumors) == "1":
        return [(n, "") for n in normals]
    paired = list(product(normals, tumors))
    if config.get("only_paired_variant_calling", True) or config.get(""):
        return paired
    return paired + [(n, "") for n in normals]


@define
class SarekTools:
    tools: list[str]
    normalize: bool


@define
class PatientOutputs:
    name: str
    tumor: str
    normal: str
    config: str
    ot: OUTPUT_TYPE = field(
        init=False,
        factory=Factory(
            lambda x: "paired"
            if x.tumor and x.normal
            else "germline"
            if x.normal
            else "tumor_only"
        ),
    )
    outdir: Path = field(
        init=False, factory=Factory(lambda x: Path(x.config["outdir"]), takes_self=True)
    )

    @property
    def sample(self) -> str:
        if self.ot == "paired":
            return f"{self.tumor}_vs_{self.normal}"
        elif self.ot == "germline":
            return self.normal
        return self.tumor

    @property
    def vep_suffix(self) -> str:
        return (
            "_VEP.ann.vcf.gz"
            if "snpeff" not in self.config["tools"]
            else "_snpEff_VEP.ann.vcf.gz"
        )

    # def bam(self)

    def vep(self, caller: str) -> str | None:
        parent = self.outdir / "annotation" / caller / self.sample
        file = parent / f"{self.sample}.{caller}_{self.vep_suffix}"
        if file.exists():
            return str(file)

    @property
    def mutect2_unfiltered(self) -> str | None:
        if self.ot == "paired":
            return None
        file = (
            self.outdir
            / "variant_calling"
            / "mutect2"
            / self.sample
            / f"{self.sample}.mutect2.vcf.gz"
        )
        return str(file)

    def all(self) -> dict:
        return {"mutect2_unfiltered": self.mutect2_unfiltered}


def rename_caller(caller: str, vcf: str) -> str:
    if caller not in {"manta", "strelka"}:
        return caller
    if caller == "manta" and "manta.diploid" in vcf:
        return "manta_diploid"
    if caller == "manta":
        return "manta_somatic"
    if caller == "strelka" and "somatic_indels" in vcf:
        return "strelka_indels"
    if caller == "strelka":
        return "strelka"


def make_sarek_manifest(config: Path, caller: str = "consensus") -> None:
    """
    Construct sarek output manifest from the yaml configuration used
    to run sarek.

    [2026-11-05 Thu] Last updated for sarek v3.11
    """
    import polars as pl

    with open(config, "r") as f:
        cfg: dict = yaml.load(f)
    tools: list[str] = cfg["tools"].split(",")
    outdir: Path = Path(cfg["outdir"])
    annotations = outdir / "annotations"

    assert "input" in cfg, "Sample manifest must be in configuration"
    assert (
        "save_mapped" in cfg and cfg["save_mapped"]
    ), "Pipeline must produce mapping files"
    assert "outdir" in cfg, "`outdir` must be in configuration"
    assert "tools" in cfg, "`tools` must be in configuration"
    assert "vep" in tools, "Pipeline must be run with vep"
    assert "mutect2" in tools, "Pipeline must be run with mutect2"
    assert caller in tools, "`caller` chosen for annotations must be in tool list"
    assert annotations.exists(), "`annotations` directory must exist"

    cfg["tools"] = tools

    T: SarekTools = SarekTools(tools=tools, normalize=cfg.get("normalize_vcfs", False))
    sample_sheet: pl.DataFrame = (
        pl.read_csv(cfg["input"], schema_overrides={"status", str})
        .with_columns(pl.struct(["status", "sample"]).alias("status_sample"))
        .group_by("patient")
        .agg(pl.col("sample_status"))
    )
    variant_called = pl.read_csv(outdir / "csv" / "variantcalled.csv")
    renamed_callers = [
        rename_caller(vc, file)
        for vc, file in zip(variant_called["variantcaller"], variant_called["vcf"])
    ]
    caller_map: dict[tuple[str, str], dict] = (
        variant_called.with_columns(pl.Series(renamed_callers).alias("variantcaller"))
        .pivot(on="variantcaller", values="vcf")
        .rows_by_key(["patient", "sample"], named=True, unique=True)
    )

    # BAM files default to tumor bam
    common = {"patient": [], "vep": [], "sample": []}

    tmp_germline = copy.deepcopy(common)
    tmp_tumor_only = copy.deepcopy(common)
    tmp_paired = copy.deepcopy(common)
    tmp_paired.update(
        {
            "bam": [],
            "bam_index": [],
            "fastq_1": [],
            "fastq_2": [],
            "mutect2_unfiltered": [],
        }
    )
    for row in sample_sheet.iter_rows(named=True):
        patient = row["patient"]
        for normal, tumor in group_samples(row["status_sample"]):
            obj = PatientOutputs(name=patient, normal=normal, tumor=tumor, config=cfg)
            vcfs: dict = caller_map[(patient, obj.sample)]
            if obj.ot == "tumor_only":
                pass
    if cfg.get("joint_germline", False) and cfg.get(
        "only_paired_variant_calling", False
    ):
        tmp_germline["patient"].append(None)
        tmp_germline["patient"].append(None)
        tmp_germline["patient"].append(None)

    if "merge" in tools:
        pass
