#!/usr/bin/env python

from SigProfilerAssignment import Analyzer as Analyze


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_folder")
    parser.add_argument("-o", "--output_folder")
    parser.add_argument(
        "-t",
        "--input_type",
        default="vcf",
        help="type of input files",
        action="store",
    )
    parser.add_argument(
        "-g",
        "--genome_build",
        default="GRCh38",
        help="Reference genome build",
        action="store",
    )
    parser.add_argument(
        "-e",
        "--exome",
        default=False,
        help="Whether or not to use exome-renormalized COSMIC signatures",
        action="store_true",
    )
    parser.add_argument(
        "-x",
        "--exclude_file",
        default="",
        help="Path to file containing cosmic signatures to exclude, each one on a newline",
        action="store",
    )
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == "__main__":
    args = parse_args()
    to_exclude: list | None = None
    if args["exclude_file"]:
        with open(args["exclude_file"], "r") as r:
            to_exclude = [s.strip() for s in r.read().splitlines()]
    Analyze.cosmic_fit(
        args["input_folder"],
        args["output_folder"],
        input_type=args["input_type"],
        exome=args["exome"],
        genome_build=args["genome_build"],
        exclude_signature_subgroups=to_exclude,
    )
