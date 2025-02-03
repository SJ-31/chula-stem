#!/usr/bin/env ipython

import re
from pathlib import Path
from typing import Callable

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from pyhere import here


def read_existing[T](
    filename: Path, expr: Callable[[Path], T], read_fn=ad.read_h5ad
) -> T:
    if filename.exists():
        return read_fn(filename)
    else:
        return expr(filename)


def match_multiple(x, patterns):
    return np.any([re.search(h, x) for h in patterns])


remote = here("analyses", "data_all")
outdir = here("analyses", "output", "corneal")

colon_dir = here(remote, "public_data", "GEO_GSE116222-human_colon")
fibro_dir = here(remote, "public_data", "GEO_GSE195452-human_skin_fibroblast")
corneal_dir = here(remote, "public_data", "GEO_GSE227942-human_corneal_endothelial")

fibro_file: Path = here(outdir, "fibro.h5ad")
corneal_file: Path = here(outdir, "corneal.h5ad")
colon_file: Path = here(outdir, "colon.h5ad")


def get_fibroblast(f):
    adata: ad.AnnData = sc.read_10x_mtx(path=here(fibro_dir, "Healthy"))
    adata.obs["source"] = "fibroblast"
    adata.obs["sample"] = "fibroblast_unknown"
    adata.write_h5ad(f)
    return adata


def get_corneal(f):
    def read(p):
        adata = sc.read_10x_mtx(path=here(corneal_dir), prefix=p)
        adata.obs["sample"] = re.sub("_.*", "", p)
        return adata

    prefixes = ("GSM7111220_g016_", "GSM7111221_g017_", "GSM7111223_g019_")
    corneal: ad.AnnData = ad.concat(
        [read(p) for p in prefixes], merge="same", index_unique="_"
    )
    corneal.obs["source"] = "corneal"
    corneal.write_h5ad(f)
    return corneal


def get_colon(f):
    healthy = ("-A1", "-C1", "-B1")
    colon = ad.AnnData(
        X=np.transpose(
            pd.read_csv(here(colon_dir, "GSE116222_Expression_matrix.txt"), sep="\t")
        )
    )
    wanted = [h for h in colon.obs.index if match_multiple(h, healthy)]
    colon = colon[wanted, :]
    colon.obs["source"] = "colon"
    colon.obs["sample"] = colon.obs.index.str.replace(".*-", "", regex=True)
    colon.write_h5ad(f)
    return colon


combined_file = here(outdir, "combined.h5ad")


def get_combined(f):
    fibro = read_existing(fibro_file, get_fibroblast)
    colon = read_existing(colon_file, get_colon)
    corneal = read_existing(corneal_file, get_corneal)
    combined = ad.concat([fibro, corneal, colon], join="outer")
    combined.write_h5ad(f)
    return combined


adata = read_existing(combined_file, get_combined)
