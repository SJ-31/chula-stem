#!/usr/bin/env python

from collections.abc import Sequence
from pathlib import Path
from typing import override

"""
Make empty files for the "wes_tumor_only" nextflow routine, allowing
 the pipeline to be started from an arbitrary entry point
"""


class SubjectResults:
    def __init__(
        self,
        subject: str,
        outdir: Path,
        sample_type: str = "tumor",
        up_to: tuple | None = None,
        make_fastq: bool = True,
    ) -> None:
        self.subject: str = subject
        self.stype: str = sample_type
        self.outdir: Path = outdir.joinpath(self.subject)
        if not self.outdir.exists():
            self.outdir.mkdir()
        if make_fastq:
            fastq = self.outdir.joinpath("fastq")
            fastq.mkdir(exist_ok=True)
            for i in ("1", "2"):
                fastq.joinpath(f"{self.subject}_{i}.fastq.gz").write_text("EMPTY")
        self.up_to: tuple | None = (
            up_to  # Tuple indicating to only fill results up to this file e.g.
            # (4, "-recal.bam") to stop making empty files
            # once 4-{subject}-recal.bam is made
        )
        calls = (
            self._fill_preprocess,
            self._fill_variant_calling,
            self._fill_annotations,
            self._fill_metrics,
        )
        self.exited: bool = False
        for call in calls:
            if call() is None:
                self.exited = True
                break

    def _mkdir(self, name: str, outdir: Path | None = None) -> Path:
        if outdir is None:
            outdir = self.outdir
        dir = outdir.joinpath(name)
        if not dir.exists():
            dir.mkdir()
        return dir

    def _fill(self, dir: str, spec: tuple) -> Path | None:
        """
        Helper function for filling a directory with empty results files

        Parameters
        ----------
        dir : directory to fill. Will be created if doesn't exist
        spec : tuple of elements (module_number, tuple of files for fill_result_folder)
        """
        path = self._mkdir(dir)
        for module_number, files in spec:
            created_files = self.fill_result_folder(
                module_number, path=path, subject=self.subject, contents=files
            )
            if self.up_to is not None and self.up_to in created_files:
                print("Finished filling up to this point")
                return None
        return path

    def _fill_preprocess(self) -> Path | None:
        dir = self._fill(
            "tumor",
            (
                (
                    1,
                    (
                        (f"_{self.stype}-fastp.json", False),
                        ("_1.fastp.fastq.gz", False),
                        ("_2.fastp.fastq.gz", False),
                        (f"_{self.stype}-fastp.html", False),
                    ),
                ),
                (
                    2,
                    ((f"_{self.stype}.bam", False),),
                ),
                (
                    3,
                    (
                        (f"_{self.stype}-dedup.bam", False),
                        (f"_{self.stype}-dedup_metrics.txt", False),
                    ),
                ),
                (
                    4,
                    (
                        (f"_{self.stype}-AnalyzeCovariates.pdf", False),
                        (f"_{self.stype}-recal.bam.bai", False),
                        (f"_{self.stype}-recal.bam_1.empty", False),
                        (f"_{self.stype}-recal_1.empty", False),
                        (f"_{self.stype}-recal.bam", False),
                    ),
                ),
            ),
        )
        if dir is not None:
            self._mkdir("4-recalibration_tables", outdir=dir)
        return dir

    def _fill_variant_calling(self) -> Path | None:
        dir = self._fill(
            "variant_calling",
            (
                (
                    5,
                    (
                        ("-candidateSmallIndels_Manta.vcf.gz", False),
                        ("-candidateSV_Manta.vcf.gz", False),
                        ("-ClairS-TO.vcf.gz", False),
                        ("-ClassifyCNV_format_cnvkit.bed", False),
                        ("-Cnvkit", True),
                        ("-Gridss_all.vcf.gz", False),
                        ("-Gridss_confident.vcf.gz", False),
                        ("-MantaOut", True),
                        ("-Msisensor", True),
                        ("-Octopus_all.vcf.gz", False),
                        ("-tumorSV_Manta.vcf.gz", False),
                    ),
                ),
            ),
        )
        if dir is not None:
            self._mkdir("5-Mutect2", dir)
        return dir

    def _fill_metrics(self) -> Path | None:
        return self._fill(
            "metrics",
            (
                (
                    8,
                    (
                        (f"_{self.stype}-Mosdepth.mosdepth.global.dist.txt", False),
                        (f"_{self.stype}-Mosdepth.mosdepth.region.dist.txt", False),
                        (f"_{self.stype}-Mosdepth.mosdepth.summary.txt", False),
                        (f"_{self.stype}-Mosdepth.per-base.bed.gz", False),
                        (f"_{self.stype}-Mosdepth.per-base.bed.gz.csi", False),
                        (f"_{self.stype}-Mosdepth.regions.bed.gz", False),
                        (f"_{self.stype}-Mosdepth.regions.bed.gz.csi", False),
                        (f"_{self.stype}-Picard_alignment_metrics.txt", False),
                        (f"_{self.stype}-Picard_hs_metrics.txt", False),
                    ),
                ),
            ),
        )

    def _fill_annotations(self) -> Path | None:
        return self._fill(
            "annotations",
            (
                (6, (("6-P12-Small_std.vcf.gz", False),)),
                (
                    7,
                    (
                        ("-ClassifyCNV_intermediate", True),
                        ("-ClassifyCNV.tsv", False),
                        ("-SigProfilerAssignment", True),
                        ("-Small_high_conf.vcf.gz", False),
                        ("-SV_high_conf.vcf.gz", False),
                        ("-VEP_small.html", False),
                        ("-VEP_small.tsv", False),
                        ("-VEP_small.vcf.gz", False),
                        ("-VEP_SV.html", False),
                        ("-VEP_SV.tsv", False),
                        ("-VEP_SV.vcf.gz", False),
                    ),
                ),
                (
                    8,
                    (
                        ("-Bcftools_stats_small", True),
                        ("-Bcftools_stats_SV", True),
                        ("-CR_CNV.tsv", False),
                        ("-CR_MSI.tsv", False),
                        ("-Small_std", True),
                        ("-SV_all", True),
                        ("-VEP_small.tsv", False),
                        ("-VEP_SV.tsv", False),
                    ),
                ),
            ),
        )

    def fill_result_folder(
        self, module_number: int | None, path: Path, subject: str, contents: Sequence
    ) -> set:
        def touch_empty(suffix: str, dir: bool = False) -> tuple:
            if module_number is not None:
                name = f"{module_number}-{subject}{suffix}"
            else:
                name = f"{subject}{suffix}"
            file: Path = path.joinpath(name)
            if file.exists():
                print(f"File {file} exists already, skipping...")
            elif dir:
                file.mkdir()
                file.joinpath("EMPTY.txt").touch()
            else:
                file.write_text("EMPTY")
            return module_number, suffix

        filled = {touch_empty(suffix, is_dir) for suffix, is_dir in contents}
        return filled


class SubjectLog(SubjectResults):
    def __init__(
        self,
        subject: str,
        outdir: Path,
        sample_type: str = "tumor",
    ) -> None:
        super().__init__(subject, outdir, sample_type, up_to=None, make_fastq=False)

    @override
    def _fill(self, dir: str, spec: tuple) -> Path | None:
        outdir: Path = self._mkdir(dir)
        for to_create in spec:
            file = outdir.joinpath(f"{to_create}.log")
            if not file.exists():
                file.write_text("EMPTY")
        return outdir

    @override
    def _fill_preprocess(self) -> Path | None:
        return self._fill("tumor", ("bqsr", "bwa", "fastp", "samtools_index"))

    @override
    def _fill_annotations(self) -> Path | None:
        outdir = self._fill(
            "annotations",
            (
                "bcftools_stats",
                "classify_cnv",
                "CNV_cross_reference",
                "concat_vcf",
                "filter_qc",
                "MSI_cross_reference",
                "qc_tsv",
                "sigprofilerassignment",
            ),
        )
        if outdir is not None:
            outdir.joinpath(f"8-{self.subject}-VEP_small.tsv_qc_tsv.log").touch()
            outdir.joinpath(f"8-{self.subject}-VEP_SV.tsv_qc_tsv.log").touch()
        return outdir

    @override
    def _fill_metrics(self) -> Path | None:
        self._fill("metrics", ("mosdepth", "picard"))

    @override
    def _fill_variant_calling(self) -> Path | None:
        return self._fill(
            "variant_calling",
            (
                "clairs-to",
                "classify_cnv_format_cnvkit",
                "cnvkit",
                "filter_mutect_calls",
                "get_pileup_summaries",
                "gridss",
                "learn_read_orientation",
                "manta",
                "MSI_cross_reference",
                "msisensor",
                "mutect2",
                "octopus",
                "run_clairs_to",
            ),
        )


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--subject_name", required=False)
    parser.add_argument("-p", "--subject_file", required=False)
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("-l", "--logdir", required=True)
    parser.add_argument("-s", "--up_to_suffix", required=False)
    parser.add_argument("-m", "--up_to_module", required=False)
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__":
    args = parse_args()
    if args["subject_name"] is None and args["subject_file"] is None:
        raise ValueError("Either subject_name or subject_file must be provided!")
    if bool(args["up_to_suffix"]) ^ bool(args["up_to_module"]):
        raise ValueError("Both up_to_suffix and up_to_module must be provided!")
    if args["up_to_suffix"] and args["up_to_module"]:
        up_to = (args["up_to_module"], args["up_to_suffix"])
    else:
        up_to = None

    outdir = Path(args["outdir"])
    outdir.mkdir(exist_ok=True)
    logdir = Path(args["logdir"])
    logdir.mkdir(exist_ok=True)

    if args["subject_file"]:
        with open(args["subject_file"], "r") as f:
            subjects = f.read().strip().splitlines()
            for s in subjects:
                sub = SubjectResults(s, outdir=outdir, up_to=up_to)
                if not sub.exited:
                    SubjectLog(s, outdir=logdir)
    else:
        sub = SubjectResults(args["subject_name"], outdir=outdir, up_to=up_to)
        if not sub.exited:
            SubjectLog(args["subject_name"], outdir=logdir)
