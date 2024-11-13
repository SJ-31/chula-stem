#!/usr/bin/env ipython


import click
import vcfpy


@click.command()
@click.option("-n", "--name", required=True, help="Name of INFO tag")
@click.option(
    "-d", "--description", required=True, help="Description of INFO tag in header"
)
@click.option(
    "-b",
    "--number",
    required=True,
    help="Number of characters included with this tag",
)
@click.option(
    "-t",
    "--type",
    required=True,
    help="Tag type (Integer, Float, Flag, Character, String)",
)
@click.option(
    "-a",
    "--default",
    required=True,
    help="Default value of tag in each record",
    type=str,
)
@click.option("-i", "--input", required=True, help="VCF file to annotate")
@click.option("-o", "--output", required=True, help="Output file")
def vcf_info_add_tag(
    name: str,
    description: str,
    number: str,
    type: str,
    default: str,
    input: str,
    output: str,
):
    _vcf_info_add_tag(name, description, number, type, default, input, output)


def _vcf_info_add_tag(
    name: str,
    description: str,
    number: str,
    type: str,
    default: str,
    input: str,
    output: str,
):
    reader = vcfpy.Reader.from_path(input)
    new_line = vcfpy.OrderedDict(
        [
            ("ID", name),
            ("Number", number),
            ("Type", type),
            ("Description", description),
        ]
    )
    reader.header.add_info_line(new_line)
    writer = vcfpy.Writer.from_path(output, reader.header)
    for record in reader:
        record: vcfpy.Record
        record.INFO[name] = [default]
        writer.write_record(record)
