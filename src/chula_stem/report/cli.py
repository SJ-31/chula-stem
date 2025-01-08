#!/usr/bin/env ipython
#

import json

import click

from chula_stem.report.variant_calling_report import VariantCallingReport


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
    paths: dict = spec.get("paths")
    if not paths:
        raise ValueError("`paths` key must be specified!")
    out = output if output else paths.get("output")
    shared_keys = {"filename", "metadata", "tmpdir"}
    kwargs: dict = {k: v for k, v in paths.items() if k not in shared_keys}
    kwargs.update(spec.get("misc", {}))
    report = None
    if report_type == "variant_calling":
        report = VariantCallingReport(
            filename=out,
            metadata=spec.get("meta", {}),
            other_text=spec.get("other_text", {}),
            tmpdir=tmpdir,
            **kwargs,
        )
    if report:
        report.build()
        report.merge()
