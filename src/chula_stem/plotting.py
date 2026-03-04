#!/usr/bin/env python
import inspect
from pathlib import Path
from typing import Literal

import anndata as ad
import click
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pygenomeviz import GenomeViz
from pygenomeviz.track import FeatureSubTrack, FeatureTrack, LinkTrack


def plot_custom_legend(ax: Axes, legend: dict, **kwargs) -> None:
    lines = [Line2D([0], [0], color=c, lw=4) for c in legend.values()]
    ax.legend(lines, list(legend.keys()), **kwargs)


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


# * pyGenomeViz tweaks


def feature_plot_all(feature: FeatureTrack, fast_render: bool = True):
    """Updated version of plot_all to use the new plot texts fn"""
    feature._plot_track_label()
    feature._plot_segment_lines()
    feature._plot_segment_sep()
    feature._plot_features(fast_render)
    feature._plot_exon_features(fast_render)
    plot_texts_wlines(feature)


def plot_texts_wlines(feature: FeatureTrack) -> None:
    """Plot texts"""
    for seg in feature.segments:
        for text_kws in seg.transform_text_kws_list:
            dct = {
                k: v for k, v in text_kws.items() if k not in {"show_line", "line_kws"}
            }
            feature.ax.text(**dct)
            if text_kws.get("show_line"):
                feature.ax.axvline(dct["x"], ymin=0.5, **text_kws.get("line_kws", {}))


def plotfig(
    gv: GenomeViz,
    *,
    dpi: int = 100,
    custom_legend: dict = None,
    fast_render: bool = True,
) -> Figure:
    """Plot figure

    Parameters
    ----------
    dpi : int, optional
        DPI
    fast_render : bool, optional
        Enable fast rendering mode using PatchCollection.
        Set fast_render=True by default, and set it to False
        when used in the `savefig_html()` method.
        Fast rendering mode cannot generate tooltips for html display.

    Returns
    -------
    fig : Figure
        Plot figure result
    """
    # Check track num
    tracks = gv.get_tracks(subtrack=True)
    if len(tracks) == 0:
        raise ValueError("Failed to plot figure. No track found!!")

    with plt.style.context(gv._mpl_style):  # type: ignore
        # Setup figure & gridspece
        fig = plt.figure(figsize=gv.figsize, dpi=dpi)
        fig.tight_layout()
        height_ratios = [t.ratio for t in tracks]
        nrows = len(tracks)
        if custom_legend:
            nrows += 2
            height_ratios.extend([0.1, 0.1])  # Dirty hack for padding between legend
        gs = GridSpec(nrows=nrows, ncols=1, height_ratios=height_ratios, hspace=0.8)
        gs.update(left=0, right=1, bottom=0, top=1, hspace=0, wspace=0)

        for idx, track in enumerate(tracks):
            # Create axes & set axes to track
            ax: Axes = fig.add_subplot(gs[idx])
            track.set_ax(ax, gv._show_axis)

            if isinstance(track, FeatureTrack):
                feature_plot_all(track)
            elif isinstance(track, FeatureSubTrack):
                pass
            elif isinstance(track, LinkTrack):
                track.plot_links(fast_render)
            else:
                track_class = track.__class__.__name__
                raise NotImplementedError(f"{track_class=} is invalid track class!!")

        lowest_track_ax = tracks[-1].ax
        if gv._plot_scale_bar:
            gv._plot_scale_bar(lowest_track_ax)

        if gv._plot_axis_ticks:
            gv._plot_axis_ticks(lowest_track_ax)

        if gv._plot_colorbar:
            gv._plot_colorbar(fig)

        if custom_legend:
            ax_new = fig.add_subplot(gs[-1])
            for pos in ("left", "right", "top", "bottom"):
                ax_new.spines[pos].set_visible(False)
            ax_new.set_zorder(tracks[-1].zorder)
            # Set xlim
            ax_new.set_xlim(tracks[-1].xlim)
            ax_new.set_ylim(*tracks[-1].ylim)
            # Set facecolor to transparent
            ax_new.set_facecolor("none")
            ax_new.tick_params(
                left=False, labelleft=False, bottom=False, labelbottom=False
            )
            fig.subplots_adjust(hspace=5)
            plot_custom_legend(ax_new, custom_legend, loc="best")

    gv._setup_jupyter_inline()

    return fig


def savefig(
    gv: GenomeViz,
    savefile: str | Path,
    *,
    custom_legend: dict = None,
    dpi: int = 100,
    pad_inches: float = 0.5,
) -> None:
    """Custom version of GenomeViz savefig to use custom plotfig"""
    fig = plotfig(gv, dpi=dpi, custom_legend=custom_legend)
    fig.savefig(
        fname=str(savefile),
        dpi=dpi,
        pad_inches=pad_inches,
        bbox_inches="tight",
    )
    if gv.clear_savefig:
        fig.clear()
        plt.close(fig)


def plot_associations(
    adata: ad.AnnData,
    groupby: str,
    assoc_df: pd.DataFrame,
    n: int,
    style: Literal["heatmap", "tracksplot"] = "tracksplot",
) -> list[Figure]:
    """
    Generate tracksplots for the top n most associated genes

    Parameters
    ----------
    assoc_df : pd.DataFrame
        DataFrame output of `find_proportional_genes`
    """
    top_assoc = (
        assoc_df.reset_index(names="genes")
        .melt(id_vars="genes")
        .groupby("variable")
        .apply(
            lambda df: df.sort_values("value", ascending=False)
            .query("genes != variable")
            .drop("variable", axis=1)
            .head(n)
        )
    ).reset_index(names=["query", "idx"])
    queries = assoc_df.columns
    result = []
    for q in queries:
        cur_genes = top_assoc.loc[top_assoc["query"] == q, :]
        if style == "tracksplot":
            vars = [q] + list(cur_genes["genes"])
            axes = sc.pl.tracksplot(
                adata,
                var_names=vars,
                show=False,
                groupby=groupby,
            )
            for ax, lab in zip(
                axes["track_axes"], ["Query"] + list(cur_genes["value"])
            ):
                ax: Axes
                old_lab = ax.get_ylabel()
                if not isinstance(lab, str):
                    lab = round(lab, 2)
                ax.set_ylabel(f"{old_lab} ({lab})")
        else:
            mapping = dict(
                Query=q,
                **{
                    str(round(v, 2)): g
                    for g, v in zip(cur_genes["genes"], cur_genes["value"])
                },
            )
            axes = sc.pl.heatmap(adata, var_names=mapping, show=False, groupby=groupby)
        result.append(axes["groupby_ax"].get_figure())
    return result
