#!/usr/bin/env ipython

from collections.abc import Callable
from pathlib import Path

import anndata as ad
import chula_stem.sc_rnaseq as sc_utils
import numpy as np
import pandas as pd
import scanpy as sc
from loguru import logger

import functions as fn

try:
    from snakemake.script import snakemake as smk
except ImportError:
    smk = type("snakemake", (), {"rule": None, "config": {}, "log": [0]})

RNG: int = smk.config["rng"]

if len(smk.log) == 1:
    logger.add(smk.log[0])

# * Functions


def call_dr(
    x: np.ndarray, method: str, kws: dict
) -> Callable[[np.ndarray], np.ndarray]:
    import torchdr as tdr

    if method == "umap":
        dr_obj = tdr.UMAP(**kws)
    elif method == "t-sne":
        dr_obj = tdr.TSNE(**kws)
    elif method == "pacmap":
        dr_obj = tdr.PACMAP(**kws)
    else:
        raise ValueError(f"DR method {method} not supported")
    return dr_obj.fit_transform(x)


# * Rules


def integrate_and_cluster(adata: ad.AnnData | None = None, cfg: dict = None):
    """
    1. Identify and subset to HVGs, accounting for batch
    2. Integrate data across batches
    3. Generate cell clusters with Leiden and assess
    """
    adata = adata or ad.read_h5ad(smk.input[0])
    cfg = cfg or smk.config.get(smk.rule, {})

    batch_key: str = cfg.get("batch_key") or smk.config["batch_key"]
    sc.pp.highly_variable_genes(
        adata, inplace=True, layer="lshift_normalized", batch_key=batch_key, subset=True
    )
    # Integration method
    i_kws = cfg["integration"]["kws"] or {}
    i_method = cfg["integration"]["method"]
    key = f"X_{i_method}"
    if i_method == "scVI":
        import scvi

        scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key)
        model = scvi.model.SCVI(adata=adata, **i_kws)
        model.train()
        key = "X_scVI"
        adata.obsm[key] = model.get_latent_representation()
    elif i_method == "harmony":
        key = "X_harmony"
        sc.external.pp.harmony_integrate(
            adata, key=[batch_key, "patient"], adjusted_basis=key
        )
    else:
        raise NotImplementedError()
    sc.external.pp.bbknn(adata, batch_key=batch_key)
    clst_cfg: dict = cfg["clustering"]
    if clst_cfg["method"] == "leiden":
        sc.pp.neighbors(adata, use_rep=key)
        sc_utils.sweep_clustering(
            adata,
            lambda x, res, key: sc.tl.leiden(
                x, resolution=res, key_added=key, **clst_cfg["kws"]
            ),
            values=clst_cfg["sweep"],
            prefix="leiden_res",
            distances=adata.obsp["distances"],
        )
    else:
        raise NotImplementedError()
    adata.write_h5ad(smk.output[0])


def prepare_data():
    _ = fn.prepare_data(smk.output[0], smk.config)


def cellassign():
    adata = ad.read_h5ad(smk.input[0])
    _, markers = fn.get_gs_and_cell_markers(smk.config)
    result = sc_utils.cell_assign_wrapper(
        adata,
        cell_markers=markers,
        model_path=Path(smk.params["model"]),
        **RCONFIG["kws"],
    )
    model = result["model"]
    result["pred"].to_csv(smk.output["predictions"], index=False)
    run_metrics = pd.concat([val for val in model.history.values()], axis=1)
    model.save(smk.params["model"], save_anndata=True, overwrite=True)
    run_metrics.reset_index().to_csv(smk.output["run_metrics"], index=False)


def do_dimensionality_reduction(unintegrated: bool = True):
    # TODO: add a modification for X to get the integrated key e.g. "X_scVI"
    cfg = smk.config["DR"]
    adata = ad.read_h5ad(smk.input[0])
    method = smk.params["method"]
    kws: dict = cfg["methods"][method].get("kws", {}) or {}
    if unintegrated:
        x: np.ndarray = adata.obsm["X_pca"][:, : cfg["n_pcs"]]
        kws["init"] = "normal"
    else:
        x = adata.obsm[smk.config["integration_key"]]
    hp_value: int | float = smk.params["hp_value"]
    hp_to_vary = cfg["methods"][method]["vary"][0]
    kws[hp_to_vary] = int(hp_value)
    result: np.ndarray = call_dr(x, method, kws)
    np.save(smk.output[0], result, allow_pickle=True)


# * Entry
if rule_fn := globals().get(smk.rule):
    rule_fn()
elif smk.rule.startswith("do_dimensionality_reduction"):
    do_dimensionality_reduction()
