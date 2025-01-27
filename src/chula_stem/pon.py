#!/usr/bin/env python

import gzip
import sys
from tempfile import NamedTemporaryFile

import click
import polars as pl


def is_gz(file: str) -> bool:
    with open(file, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def parse_vcf_header(file: str) -> tuple[str, list]:
    header_lines: list[str] = []

    def check_header(text, index: int) -> list:
        header_line = (
            "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"
        )
        if header_line in text:
            samples = text.replace(header_line, "").strip().split("\t")
            return samples
        elif "#" not in text and index != 0:
            raise ValueError("No sample found")
        return []

    open_fn = gzip.open if is_gz(file) else open
    lines = open_fn(file, "rb")
    for index, line in enumerate(lines):
        text = line.decode().strip()
        checked = check_header(text, index)
        header_lines.append(text)
        if checked:
            header_text = f"{'\n'.join(header_lines)}\n"
            return header_text, checked
    return "", []


def get_schema(samples: list) -> tuple[dict, list]:
    schema = {
        "CHROM": pl.String,
        "POS": pl.Int64,
        "ID": pl.String,
        "REF": pl.String,
        "ALT": pl.String,
        "QUAL": pl.String,
        "FILTER": pl.String,
        "INFO": pl.String,
        "FORMAT": pl.String,
    }
    for s in samples:
        schema[s] = pl.String
    return schema, list(schema.keys())


@click.command()
@click.option("-f", "--file", required=False, help="Input vcf file", default="")
@click.option("-o", "--output", required=False, help="Output file", default="")
@click.option(
    "-m",
    "--minimum",
    required=False,
    help="Minimum number of samples variant must be found in",
    default=2,
)
def filter_minimum_samples(file: str = "", minimum: int = 2, output: str = "") -> None:
    from_stdin: bool = not file
    if from_stdin:
        tmp = NamedTemporaryFile("w+b", delete_on_close=True, dir=".")
        tmp.write(sys.stdin.read().encode())
        input: str = tmp.name
    else:
        input = file
    header, samples = parse_vcf_header(input)
    schema, cols = get_schema(samples)
    df = pl.read_csv(
        input,
        separator="\t",
        comment_prefix="#",
        schema=schema,
        new_columns=cols,
        null_values=["."],
    )
    # Count the number of samples that the variant occurs in
    counts: pl.Series = (
        df.select(samples)
        .with_columns(
            pl.col(samples).str.replace_all("[\\./:]", "").str.len_chars() > 0
        )
        .cast(pl.Int64)
        .sum_horizontal()
    )
    df = df.filter(counts >= minimum)
    if not output:
        sys.stdout.write(header)
        df.write_csv(sys.stdout, separator="\t", null_value=".", include_header=False)
    else:
        with open(output, "a") as f:
            f.write(header)
            df.write_csv(f, separator="\t", null_value=".", include_header=False)
    if from_stdin:
        tmp.close()
