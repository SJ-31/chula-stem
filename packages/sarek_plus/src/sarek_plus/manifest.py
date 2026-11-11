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


def rename_caller(caller: str, vcf: str) -> str:
    if caller not in {"manta", "strelka", "tiddit"}:
        return caller
    if caller == "tiddit" and "germline.vcf.gz" in vcf:
        return "tiddit_germline"
    if caller == "tiddit":
        return caller
    if caller == "manta" and "manta.diploid" in vcf:
        return "manta_diploid"
    if caller == "manta":
        return "manta_somatic"
    if caller == "strelka" and "somatic_indels" in vcf:
        return "strelka_indels"
    if caller == "strelka":
        return "strelka"


def get_fq(
    sample, patient, outdir: Path, config: dict, fastq_map: dict
) -> tuple[str, str | None]:
    """
    If fastq files were split, the first element is the directory containing
    the fastp files and the second is None.
    Otherwise, return a tuple of fastq1, fastq2
    """
    if "save_trimmed" in config:
        dir: Path = outdir / "preprocessing" / "fastp" / sample
        if len(list(dir.iterdir())) > 2:
            return str(dir), None
        return str(next(dir.glob("*R1.fastp.fastq.gz"))), str(
            next(dir.glob("*R2.fastp.fastq.gz"))
        )
    fq1 = fastq_map.get((patient, sample))["fastq_1"]
    fq2 = fastq_map.get((patient, sample))["fastq_2"]
    return fq1, fq2


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


def _make_sarek_manifest(config: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    """
    Construct sarek output manifest from the yaml configuration used
    to run sarek.

    [2026-11-05 Thu] Last updated for sarek v3.11

    Returns
    -------
    Tuple of dataframes: (tool_outs, preprocessed)
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
    assert (
        "save_output_as_bam" in cfg and cfg["save_output_as_bam"]
    ), "Pipeline must save mapping files as BAM"
    assert "outdir" in cfg, "`outdir` must be in configuration"
    assert "tools" in cfg, "`tools` must be in configuration"
    assert "vep" in tools, "Pipeline must be run with vep"
    assert "mutect2" in tools, "Pipeline must be run with mutect2"
    assert not cfg.get("trim_fastq", False) or (
        cfg["trim_fastq"] and cfg.get("save_trimmed", False)
    ), "If trimming is enabled, must save trimmed output"
    assert annotation.exists(), "`annotations` directory must exist"

    cfg["tools"] = tools

    sample_sheet = pl.read_csv(cfg["input"], schema_overrides={"status": str})
    samples_grouped: pl.DataFrame = (
        sample_sheet.with_columns(
            pl.struct(["status", "sample"]).alias("status_sample")
        )
        .group_by("patient")
        .agg(pl.col("status_sample"))
    )
    fq_map: dict = sample_sheet.select(
        ["patient", "sample", "fastq_1", "fastq_2"]
    ).rows_by_key(["patient", "sample"], named=True, unique=True)

    variant_called = pl.read_csv(outdir / "csv" / "variantcalled.csv")
    if cfg.get("normalize_vcfs"):
        # NOTE: [2026-11-11 Wed] For some reason the "variantcalled" file
        # includes the normalized versions of the vcf file and the
        # unnormalized ones, but inconsistently for each caller
        variant_called = (
            variant_called.group_by("patient", "sample", "variantcaller")
            .agg("vcf")
            .with_columns(
                pl.when(
                    pl.col("vcf")
                    .list.eval(pl.element().str.ends_with(".norm.sorted.vcf.gz"))
                    .list.any()
                )
                .then(
                    pl.col("vcf").list.filter(
                        pl.element().str.ends_with(".norm.sorted.vcf.gz")
                    )
                )
                .otherwise(pl.col("vcf"))
            )
            .explode("vcf")
        )
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
    tool_outs, preprocessed = defaultdict(list), defaultdict(list)
    for row in samples_grouped.iter_rows(named=True):
        patient = row["patient"]
        added_preprocessed: bool = False
        for normal, tumor in group_samples(row["status_sample"], cfg):
            obj = PatientOutputs(name=patient, normal=normal, tumor=tumor, config=cfg)

            tool_outs["output_type"].append(obj.ot)
            tool_outs["normal_sample"].append(obj.normal)
            tool_outs["tumor_sample"].append(obj.tumor)
            tool_outs["id"].append(patient)
            tool_outs["mutect2_unfiltered"].append(obj.mutect2_unfiltered)
            vcfs: dict = caller_map[(patient, obj.sample)]
            for caller, file in vcfs.items():
                if (
                    caller == "haplotypecaller"
                    and joint_germline
                    and obj.ot == "germline"
                ):
                    file = joint_germline
                annotated: str | None = obj.vep(file)
                tool_outs[caller].append(file)
                tool_outs[f"vep_{caller}"].append(annotated)

            if not added_preprocessed:
                tumor_recal, normal_recal = (
                    recal.get((patient, obj.tumor), {}),
                    recal.get((patient, obj.normal), {}),
                )
                preprocessed["id"].extend([patient, patient])
                preprocessed["status"].extend(["1", "0"])
                preprocessed["recal"].extend(
                    [tumor_recal.get("bam"), normal_recal.get("bam")]
                )
                preprocessed["recal_index"].extend(
                    [tumor_recal.get("bai"), normal_recal.get("bai")]
                )
                n_fastq_1, n_fastq_2 = get_fq(obj.normal, patient, outdir, cfg, fq_map)
                t_fastq_1, t_fastq_2 = get_fq(obj.tumor, patient, outdir, cfg, fq_map)
                preprocessed["fastq_1"].extend([t_fastq_1, n_fastq_1])
                if not n_fastq_2:
                    preprocessed["fastq_2"].extend([None, None])
                    preprocessed["split_fastq"].extend([True, True])
                else:
                    preprocessed["fastq_2"].extend([t_fastq_2, n_fastq_2])
                    preprocessed["split_fastq"].extend([False, False])
                added_preprocessed = True

    return pl.DataFrame(tool_outs), pl.DataFrame(preprocessed)


@click.command()
@click.option(
    "-c",
    "--config",
    required=True,
    help="YAML file used in original sarek run",
    default=None,
)
@click.option(
    "-p",
    "--preprocessed",
    required=True,
    help="CSV file for preprocessed data",
    default=None,
)
@click.option(
    "-t",
    "--tool_outs",
    required=True,
    help="CSV file for tumor-only output",
    default=None,
)
def make_sarek_manifest(config, preprocessed, tool_outs) -> None:
    t, p = _make_sarek_manifest(config=config)
    t = pl.DataFrame({"id": []}) if t.is_empty() else t
    p = pl.DataFrame({"id": []}) if p.is_empty() else p
    t.write_csv(tool_outs)
    p.write_csv(preprocessed)


def test_manifest():
    from pyhere import here

    test_config = here("nextflow", "input", "manifests", "sarek_test_pdac_1.yaml")
    tool_outs, pp = _make_sarek_manifest(test_config)
    tool_outs.write_csv(here("tests", "sarek_plus_out") / "tool_outs_manifest.csv")
    pp.write_csv(here("tests", "sarek_plus_out") / "preprocessed_manifest.csv")
