#!/usr/bin/env ipython

from collections import defaultdict
from itertools import product
from pathlib import Path
from typing import Literal, TypeAlias

import click
import polars as pl
import yaml
from attrs import Factory, define, field

OUTPUT_TYPE: TypeAlias = Literal["tumor_only", "paired", "germline"]


@define
class PatientOutputs:
    name: str
    tumor: str
    normal: str
    config: str
    ot: OUTPUT_TYPE = field(
        init=False,
        default=Factory(
            lambda x: "paired"
            if x.tumor and x.normal
            else "germline"
            if x.normal
            else "tumor_only",
            takes_self=True,
        ),
    )
    outdir: Path = field(
        init=False, default=Factory(lambda x: Path(x.config["outdir"]), takes_self=True)
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

    def vep(self, variant_vcf: str) -> str | None:
        if variant_vcf:
            vcf_path = Path(variant_vcf)
            name = vcf_path.stem.removesuffix(".vcf")
            caller = vcf_path.parent.parent.stem
            dir = self.outdir / "annotation" / caller / self.sample
            file = dir / f"{name}{self.vep_suffix}"
            if file.exists():
                return str(file)

    @property
    def mutect2_unfiltered(self) -> str | None:
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


def drop_null_cols(df: pl.DataFrame) -> pl.DataFrame:
    nulls = [col for col, series in zip(df.columns, df) if series.is_null().all()]
    return df.drop(nulls)


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
    if config.get("only_paired_variant_calling", True):
        return paired
    return paired + [(n, "") for n in normals]


def _make_sarek_manifest(
    config: Path,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """
    Construct sarek output manifest from the yaml configuration used
    to run sarek.

    [2026-11-05 Thu] Last updated for sarek v3.11

    Returns
    -------
    Tuple of dataframes: (germline, tumor_only, paired)
    """

    with open(config, "r") as f:
        cfg: dict = yaml.full_load(f)
    tools: list[str] = cfg["tools"].split(",")
    outdir: Path = Path(cfg["outdir"])
    annotation = outdir / "annotation"

    assert "input" in cfg, "Sample manifest must be in configuration"
    assert (
        "save_mapped" in cfg and cfg["save_mapped"]
    ), "Pipeline must produce mapping files"
    assert "outdir" in cfg, "`outdir` must be in configuration"
    assert "tools" in cfg, "`tools` must be in configuration"
    assert "vep" in tools, "Pipeline must be run with vep"
    assert "mutect2" in tools, "Pipeline must be run with mutect2"
    assert annotation.exists(), "`annotations` directory must exist"

    cfg["tools"] = tools

    sample_sheet: pl.DataFrame = (
        pl.read_csv(cfg["input"], schema_overrides={"status": str})
        .with_columns(pl.struct(["status", "sample"]).alias("status_sample"))
        .group_by("patient")
        .agg(pl.col("status_sample"))
    )
    variant_called = pl.read_csv(outdir / "csv" / "variantcalled.csv")
    recal: dict[tuple(str, str), dict] = pl.read_csv(
        outdir / "csv" / "recalibrated.csv"
    ).rows_by_key(["patient", "sample"], unique=True, named=True)

    renamed_callers = [
        rename_caller(vc, file)
        for vc, file in zip(variant_called["variantcaller"], variant_called["vcf"])
    ]

    caller_map: dict[tuple[str, str], dict] = (
        variant_called.with_columns(pl.Series(renamed_callers).alias("variantcaller"))
        .pivot(on="variantcaller", values="vcf")
        .rows_by_key(["patient", "sample"], named=True, unique=True)
    )
    joint_germline: str | None = caller_map.pop(
        ("all_samples", "joint_variant_calling"), {}
    ).get("haplotypecaller")

    # BAM files default to tumor bam
    germline, t_only, paired = defaultdict(list), defaultdict(list), defaultdict(list)
    for row in sample_sheet.iter_rows(named=True):
        patient = row["patient"]
        for normal, tumor in group_samples(row["status_sample"], cfg):
            obj = PatientOutputs(name=patient, normal=normal, tumor=tumor, config=cfg)
            tumor_recal, normal_recal = (
                recal.get((patient, obj.tumor), {}),
                recal.get((patient, obj.normal), {}),
            )
            if obj.ot == "paired":
                to_update = paired
                to_update["mutect2_unfiltered"].append(obj.mutect2_unfiltered)
                to_update["recal_tumor"].append(tumor_recal.get("cram"))
                to_update["recal_tumor_index"].append(tumor_recal.get("crai"))
                to_update["recal_normal"].append(normal_recal.get("cram"))
                to_update["recal_normal_index"].append(normal_recal.get("crai"))
            elif obj.ot == "tumor_only":
                to_update = t_only
                to_update["recal"].append(tumor_recal.get("cram"))
                to_update["recal_index"].append(tumor_recal.get("crai"))
            else:
                to_update = germline
                to_update["recal"].append(normal_recal.get("cram"))
                to_update["recal_index"].append(normal_recal.get("crai"))
            to_update["id"].append(patient)
            vcfs: dict = caller_map[(patient, obj.sample)]
            for caller, file in vcfs.items():
                if caller == "haplotypecaller" and joint_germline:
                    file = joint_germline
                annotated: str | None = obj.vep(file)
                to_update[caller].append(file)
                to_update[f"vep_{caller}"].append(annotated)

    return (
        pl.DataFrame(germline).pipe(drop_null_cols),
        pl.DataFrame(t_only).pipe(drop_null_cols),
        pl.DataFrame(paired).pipe(drop_null_cols),
    )


@click.command()
@click.option(
    "-c",
    "--config",
    required=True,
    help="YAML file used in original sarek run",
    default=None,
)
@click.option(
    "-g",
    "--germline",
    required=True,
    help="CSV file for germline output",
    default=None,
)
@click.option(
    "-t",
    "--tumor_only",
    required=True,
    help="CSV file for tumor-only output",
    default=None,
)
@click.option(
    "-p", "--paired", required=True, help="CSV file for paired output", default=None
)
def make_sarek_manifest(config, germline, tumor_only, paired) -> None:
    g, t, p = _make_sarek_manifest(config=config)
    g.write_csv(germline)
    t.write_csv(tumor_only)
    p.write_csv(paired)


def test_manifest():
    from pyhere import here

    test_config = here("nextflow", "input", "manifests", "sarek_test_pdac_1.yaml")
    germline, t_only, paired = _make_sarek_manifest(test_config)
    germline.write_csv(here("tests", "sarek_plus_out") / "germline_manifest.csv")
    # t_only.write_csv(here("tests", "sarek_plus_out") / "t_only_manifest.csv")
    paired.write_csv(here("tests", "sarek_plus_out") / "paired_manifest.csv")
