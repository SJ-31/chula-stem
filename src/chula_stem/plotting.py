#!/usr/bin/env python
import rpy2
import click
from rpy2.robjects.packages import STAP
from rpy2.robjects import NULL, ListVector
from pathlib import Path
import inspect


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


def none2null(item):
    if item == None:
        return NULL
    if isinstance(item, dict):
        changed = dict(
            list(map(lambda x: (x[0], NULL) if x[1] == None else x, item.items()))
        )
        return changed
    if isinstance(item, tuple):
        change_fn = tuple
    else:
        change_fn = list
    return change_fn(map(lambda x: NULL if x == None else x, item))


def plot_cnvkit(cnr: str, cns: str, chr: str, sizing: dict, output: str = "cnvkit.png"):
    source: str = get_rscripts().joinpath("plotting.R").read_text()
    chr = NULL if not chr else chr
    conv = ListVector(none2null(sizing))
    plot_lib: STAP = STAP(source, "plotting")
    plot_lib.plot_cnvkit(cns, cnr, chr.split(","), conv, output)


def main(args, sizing):
    if args["name"] == "cnvkit":
        plot_cnvkit(
            cnr=args["cnr"],
            cns=args["cns"],
            chr=args["loci"],
            output=args["output"],
            sizing=sizing,
        )


def get_rscripts() -> Path:
    src = Path(inspect.getfile(inspect.currentframe())).parent.parent
    return src.joinpath("R")


if __name__ == "__main__":
    args, sizing = parse_args()
    main(args, sizing)
