#!/usr/bin/env python
import inspect
from pathlib import Path

import click


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="Plot name", dest="name")
    parser.add_argument("-o", "--output", help="Name of output plot file")
    parser.add_argument("-d", "--dpi", help="Plot DPI", required=False, default=300)
    parser.add_argument(
        "-w", "--width", help="Width of plot", required=False, default=None
    )
    # In matplotlib, must set on the figure e.g. figure.set_size_inches()
    parser.add_argument(
        "-u",
        "--units",
        help="Units of width and height",
        required=False,
        default="in",
    )
    parser.add_argument(
        "-h", "--height", help="Height of plot", required=False, default=None
    )
    cnvkit = subparsers.add_parser(
        "cnvkit", help="Plot copy number ratios/integer values from cnvkit data"
    )
    cnvkit.add_argument("-c", "--cns", help="Segmentation file")
    cnvkit.add_argument("-r", "--cnr", help="Copy number ratio file")
    cnvkit.add_argument(
        "-l",
        "--loci",
        help="Loci to plot, a comma-delimited list e.g. chr1,chr2 or with ranges chr1:100-100000",
        required=False,
        default=None,
    )
    args = vars(parser.parse_args())  # convert to dict
    sizing = {
        "width": args["width"],
        "height": args["height"],
        "dpi": args["dpi"],
        "units": args["units"],
    }
    return args, sizing


def get_rscripts() -> Path:
    src = Path(inspect.getfile(inspect.currentframe())).parent.parent
    return src.joinpath("R")


if __name__ == "__main__":
    args, sizing = parse_args()
    plot_main(args, sizing)
