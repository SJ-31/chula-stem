#!/usr/bin/env ipython

from typing import Callable

import anndata as ad
import numpy as np
from loguru import logger

import functions as fn

try:
    from snakemake.script import snakemake as smk
except ImportError:
    smk = type("snakemake", (), {"rule": None, "config": {}, "log": [0]})

RCONFIG = smk.config(smk.rule)
RNG: int = smk.config["rng"]

logger.add(smk.log[0])

# * Functions


def call_dr(
    x: np.ndarray, method: str, kws: dict
) -> Callable[[np.ndarray], np.ndarray]:
    import torchdr as tdr

    if method == "umap":
        dr_obj = tdr.UMAP
    elif method == "t-sne":
        dr_obj = tdr.TSNE
    else:
        raise ValueError(f"DR method {method} not supported")
    return dr_obj.fit_transform(x, **kws)


# * Rules


def prepare_data():
    _ = fn.prepare_data(smk.output[0], smk.config)


def do_dimensionality_reduction(adata_path: ad.AnnData):
    cfg = smk.config["DR"]
    adata = ad.read_h5ad(adata_path)
    x: np.ndarray = adata.obsm["X_pca"][: cfg["n_pcs"]]
    hp_value: int | float = smk.params["hp_value"]
    method = smk.params["method"]
    kws: dict = cfg[method].get("kws", {}) or {}
    hp_to_vary = cfg[method]["vary"][0]
    kws[hp_to_vary] = int(hp_value)
    result: np.ndarray = call_dr(x, method, kws)
    np.save(smk.output[0], result, allow_pickle=True)


# Call function named after rule automatically
# * Entry
if not (fn := globals().get(smk.rule)):
    raise ValueError(
        f"Function for rule Symbol’s value as variable is void: {smk.rule} not defined in this file"
    )
else:
    fn()
