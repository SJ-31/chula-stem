#!/usr/bin/env ipython
import importlib.resources as res
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal

import anndata as ad
import chula_stem.r_utils as ru
import chula_stem.sc_rnaseq as sc_utils
import chula_stem.utils as cu
import decoupler as dc
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotnine as gg
import polars as pl
import pymupdf
import rpy2.robjects as ro
import scanpy as sc
import scipy.spatial.distance as spd
import yaml
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.text import Annotation, Text
from pyhere import here
from pymupdf import Document
from scanpy.plotting import DotPlot
from yte import process_yaml

matplotlib.use("QtAgg")

workdir = here("analyses", "scrnaseq_hcc")


sys.argv.append("test=True")
with open(workdir / "cellranger_config.yaml", "r") as f:
    env = process_yaml(f)


def load_file_as_module(name, location):
    from importlib import util

    spec = util.spec_from_file_location(name, location)
    module = util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


fn = load_file_as_module("functions", workdir / "functions.py")
# %%

adata = ad.read_h5ad(workdir / "test" / "filtered.h5ad")
adata = adata[:, ~adata.var["hgnc_symbol"].isnull()].copy()

df = pl.from_pandas(adata.obs)

# agg_df = adata.obs.groupby("sample").agg("first")
# agg = sc.get.aggregate(adata, "sample", "sum")

# * Enrichment analyses

clusters = ad.read_h5ad(workdir / "test" / "seurat" / "harmony_clustering.h5ad")

gcol = "leiden_res1"
tmp = adata.copy()
tmp.obs = tmp.obs.merge(clusters.obs, left_index=True, right_index=True, how="left")
tmp = tmp[~tmp.obs[gcol].isna(),]

gs_df, _ = fn.get_gs_and_cell_markers(env, (5, None), True, ("source", "target"))
gs_df = gs_df.drop_duplicates()

dc.mt.decouple(tmp, gs_df, methods=["ulm"])
scores = dc.pp.get_obsm(tmp, key="score_ulm")
ranked = dc.tl.rankby_group(scores, groupby=gcol)


df = pd.read_csv(
    "/home/shannc/Bio_SDD/chula-stem/analyses/data_all/output/HCC/SCRNASEQ/cohort/annotations/clusters-edgeR_de.csv"
)

dummy = ad.AnnData(X=np.zeros_like(adata.X.toarray()), var=adata.var)
tmp = sc_utils.find_proportional_genes_rho(dummy, ["A1CF", "A2M", "A1BG-AS1"])

axes = sc.pl.heatmap(
    adata,
    var_names={"foo": adata.var_names[:5], "bar": adata.var_names[10:15]},
    groupby="sample",
    show=False,
)
