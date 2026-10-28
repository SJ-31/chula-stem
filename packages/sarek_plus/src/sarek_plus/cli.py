#!/usr/bin/env ipython
#

import json

import click

from sarek_plus.variant_calling_report import VariantCallingReport


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
