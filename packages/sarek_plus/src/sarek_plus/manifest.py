#!/usr/bin/env ipython

from collections import defaultdict
from itertools import product
from pathlib import Path
from typing import Literal, TypeAlias

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

    def vep(self, variant_vcf: str) -> str | None:
        vcf_path = Path(variant_vcf)
        name = vcf_path.stem.removesuffix(".vcf")
        caller = vcf_path.parent.parent.stem
        dir = self.outdir / "annotation" / caller / self.sample
        file = dir / f"{name}_{self.vep_suffix}"
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


def make_sarek_manifest(
    config: Path, caller: str = "consensus"
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
    joint_germline: str | None = renamed_callers.pop(
        ("all_samples", "joint_variant_calling"), {}
    ).get("haplotypecaller")

    caller_map: dict[tuple[str, str], dict] = (
        variant_called.with_columns(pl.Series(renamed_callers).alias("variantcaller"))
        .pivot(on="variantcaller", values="vcf")
        .rows_by_key(["patient", "sample"], named=True, unique=True)
    )

    # BAM files default to tumor bam
    germline, t_only, paired = defaultdict(list), defaultdict(list), defaultdict(list)
    for row in sample_sheet.iter_rows(named=True):
        patient = row["patient"]
        for normal, tumor in group_samples(row["status_sample"]):
            obj = PatientOutputs(name=patient, normal=normal, tumor=tumor, config=cfg)
            if obj.ot == "paired":
                to_update = paired
                to_update["mutect2_unfiltered"].append(obj.mutect2_unfiltered)
            elif obj.ot == "tumor_only":
                to_update = t_only
            else:
                to_update = germline
            to_update["patient"].append(patient)
            vcfs: dict = caller_map[(patient, obj.sample)]
            for caller, file in ((c, v) for c, v in vcfs.items() if v):
                if caller == "haplotypecaller" and joint_germline:
                    file = joint_germline
                annotated: str | None = obj.vep(file)
                to_update[caller].append(file)
                to_update[f"vep_{caller}"].append(annotated)

    return pl.DataFrame(germline), pl.DataFrame(t_only), pl.DataFrame(paired)
