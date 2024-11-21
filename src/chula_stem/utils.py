#!/usr/bin/env ipython


from tempfile import TemporaryFile

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
    output: str,  # Uncompressed vcf only
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
