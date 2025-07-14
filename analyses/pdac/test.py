#!/usr/bin/env ipython

import marsilea as ma
import marsilea.layers as mal
import marsilea.plotter as mapl
import marsilea.utils as mau
import numpy as np
import pandas as pd
from oncoprinter import OncoPrint
from oncoprinter.core import GenomicData, _format_percentage
from pyhere import here

outdir = here("analyses", "output", "pdac")

df = pd.read_csv(here(outdir, "replicate_fig_test.csv"))

conse = df["Consequence"].value_counts()

# [2025-06-05 Thu] TODO: main thing to do is modifying the GenomicData class
# to return custom pieces and not whatever the authors define
# If you want to make a custom oncoprint, define your own GenomicData class
# or subclass the authors'
# include more variant types and pass it to oncoprint here
# See https://github.com/Marsilea-viz/marsilea/blob/main/src/oncoprinter/core.py


class OP(ma.ClusterBoard):
    """OncoPrint plot

    The oncoprint plot is a visualization for genomics data in cancer research.
    It's first introduced by the cBioPortal project.
    See https://www.cbioportal.org/oncoprinter for more details.

    To use this class, import from oncoprinter

        >>> from oncoprinter import OncoPrint

    Parameters
    ----------
    genomic_data : pd.DataFrame
        Genomics data, each column is:
            1) Sample ID
            2) Track name
            3) Alteration

    patients_order : list, optional
        The order of samples, by default None
    tracks_order : list, optional
        The order of tracks, by default None
    pieces : dict, optional
        Custom pieces for each alteration, by default None
        See :class:`Piece <marsilea.layers.Piece>` for details
    background_color : str, optional, default: "#BEBEBE"
        The background color
    shrink : tuple, optional, default: (0.8, 0.8)
        The shrink ratio for each layer
    width, height : float, optional
        The size in inches to define the size of main canvas
    aspect : float, optional, default: 2.5
        The aspect ratio of the main canvas
    legend_kws : dict, optional
        The options for legend, by default None
        See :class:`cat_legend <legendkit.cat_legend>` for details
    name : str, optional
        The name of this OncoPrint
    add_tracks_names : str, optional, default: "left"
        The position to add tracks names
        If None, will not add tracks names
    add_samples_names : str, optional, default: "bottom"
        The position to add samples names
        If None, will not add samples names
    add_mut_perc : str, optional, default: "right"
        The position to add mutation percentage
        If None, will not add mutation percentage
    add_tracks_counts : str, optional, default: "right"
        The position to add tracks mutation counts
        If None, will not add tracks mutation counts
    add_mut_counts : str, optional, default: "top"
        The position to add mutation counts
        If None, will not add mutation counts
    add_tracks_counts_size : float, optional, default: 0.2
        The size of tracks mutation counts
    add_tracks_counts_pad : float, optional, default: 0
        The padding of tracks mutation counts
    add_mut_counts_size : float, optional, default: 0.2
        The size of mutation counts
    add_mut_counts_pad : float, optional, default: 0.1
        The padding of mutation counts

    """

    def __init__(
        self,
        genomic_data=None,
        patients_order=None,
        tracks_order=None,
        pieces=None,
        background_color="#BEBEBE",
        shrink=(0.8, 0.8),
        width=None,
        height=None,
        aspect=2.5,
        legend_kws=None,
        name=None,
        add_tracks_names="left",
        add_samples_names="bottom",
        add_mut_perc="right",
        add_tracks_counts="right",
        add_tracks_counts_size=0.2,
        add_tracks_counts_pad=0,
        add_mut_counts_size=0.2,
        add_mut_counts_pad=0.1,
    ):
        data = GenomicData(
            genomic_data,
            samples_order=patients_order,
            tracks_order=tracks_order,
            custom_pieces=pieces,
        )
        self.genomic_data = data
        width, height = mau.get_canvas_size_by_data(
            data.shape, width=width, height=height, scale=0.2, aspect=aspect
        )

        self.canvas = super().__init__(
            name=name, cluster_data=np.zeros(data.shape), width=width, height=height
        )

        legend_options = dict(title="Alterations", handleheight=aspect, handlelength=1)
        legend_kws = {} if legend_kws is None else legend_kws
        legend_options.update(legend_kws)

        layers, pieces, colors_mapper = [], [], {}
        for layer in data.get_layers_data(background_color):
            layers.append(layer.matrix)
            pieces.append(layer.piece)
            colors_mapper[layer.event] = layer.color
        mesh = mal.LayersMesh(
            layers=layers, pieces=pieces, shrink=shrink, legend_kws=legend_options
        )
        self.add_layer(mesh)

        if add_tracks_names:
            self.add_plot(add_tracks_names, mapl.Labels(data.tracks))

        # Add other statistics
        track_mut_rate = data.get_track_mutation_rate()
        # Convert to percentage string
        if add_mut_perc:
            rates = [_format_percentage(t) for t in track_mut_rate]
            self.add_plot(add_mut_perc, mapl.Labels(rates))
        if add_samples_names:
            self.add_plot(add_samples_names, mapl.Labels(data.samples))

        if add_tracks_counts:
            track_counter = data.get_track_mutation_types()
            colors = [colors_mapper[e] for e in track_counter.index]
            track_bar = mapl.StackBar(track_counter, colors=colors, show_value=False)
            self.add_plot(
                add_tracks_counts,
                track_bar,
                legend=False,
                size=add_tracks_counts_size,
                pad=add_tracks_counts_pad,
            )

        self.add_legends()

    clinical_plots = {
        "bar": mapl.Numbers,
        "stack_bar": mapl.StackBar,
    }


fig = OP(
    df.loc[:, ["sample", "SYMBOL", "Consequence"]],
    pieces={
        "stop_gained": mal.FrameRect(color="#7E2E84", width=2),
        "intron_variant": mal.Rect(color="#D14081"),
        "missense_variant": mal.RightTri(color="#EF798A"),
        "3_prime_UTR_variant": mal.FracRect(color="#F9F5E3", frac=(0.5, 0.5)),
        "frameshift_variant": mal.RightTri(color="#CCF5AC", right_angle="upper right"),
    },
    tracks_order=["KRAS", "TP53", "SMAD4", "CDKN2A", "KMT2C"],
)
fig.render()
fig.figure.show()

onco = ma.load_data("oncoprint")
