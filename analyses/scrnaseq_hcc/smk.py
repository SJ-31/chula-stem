#!/usr/bin/env ipython

from collections.abc import Callable
from pathlib import Path

import anndata as ad
import chula_stem.sc_rnaseq as sc_utils
import numpy as np
import pandas as pd
import plotnine as gg
import scanpy as sc
from loguru import logger

import functions as fn

try:
    from snakemake.script import snakemake as smk
except ImportError:
    smk = type("snakemake", (), {"rule": None, "config": {}, "log": [0]})

RNG: int = smk.config["rng"]
RCONFIG: dict = smk.config.get(smk.rule) or {}

if len(smk.log) == 1:
    logger.add(smk.log[0])


# * Utility functions


def get_integration_key(integration_method: str) -> str:
    return f"X_{integration_method}"


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


def integrate_data(adata: ad.AnnData, batch_key, name: str, extras: dict, env) -> str:
    params = env["integration"][name]
    method = params.get("method", name)
    key = get_integration_key(method)
    kws = params.get("kws") or {}
    if method == "scVI":
        import scvi

        scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key)
        train_kws = kws.pop("train_kws", {})
        model = scvi.model.SCVI(adata=adata, **kws)
        model.train(**train_kws)
        adata.obsm[key] = model.get_latent_representation()
        extras["model"] = model
    elif method == "harmony":
        import harmonypy as hm

        # Re-run PCA with a different feature set
        sc.pp.pca(adata, layer=params.get("layer"))
        harmony_out = hm.run_harmony(adata.obsm["X_pca"], adata.obs, vars_use=batch_key)
        adata.obsm[key] = harmony_out.Z_corr
    else:
        raise NotImplementedError()
    return key


def select_features(
    name: str, adata: ad.AnnData, batch_key: str, layer: str, env: dict
):
    cur_cfg = env["feature_selection"][name] or {}
    method = cur_cfg.get("method", name)
    kws = cur_cfg.get("kws") or {}
    selection_fn = fn.fs_dispatch(method, adata)
    if method == "seurat":
        kws.update(
            {
                "inplace": True,
                "layer": layer,
                "batch_key": batch_key,
                "subset": True,
            }
        )
        selection_fn(**kws)
    else:
        raise NotImplementedError()
    return adata


# * Rules


def integrate(adata: ad.AnnData | None = None, cfg: dict = None):
    """
    1. Identify and subset to HVGs, accounting for batch
    2. Integrate data across batches
    3. Generate cell clusters with Leiden and assess
    """
    adata = adata or ad.read_h5ad(smk.input[0])
    cfg = cfg or RCONFIG

    batch_key: str = cfg.get("batch_key") or smk.config["batch_key"]
    fs_name = smk.params["feature_selection"]
    adata = select_features(
        name=fs_name,
        adata=adata,
        env=smk.config,
        layer="lshift_normalized",
        batch_key=batch_key,
    )

    integration_name = smk.params["integration"]
    extras = {}
    key = integrate_data(
        adata, batch_key, name=integration_name, extras=extras, env=smk.config
    )
    if "model" in extras:
        extras["model"].save(
            **fn.get_scvi_model_path(fs_name, integration_name, smk.config)
        )
    integrated = ad.AnnData(obs=adata.obs, obsm={key: adata.obsm[key]}, var=adata.var)
    del integrated.obs
    del integrated.var
    integrated.write_h5ad(smk.output[0])


def cluster_cells(cfg: dict = None):
    adata = ad.read_h5ad(smk.input[0])
    integration_layer = ad.read_h5ad(smk.input[1])
    adata.obsm.update(integration_layer.obsm)
    cfg = cfg or RCONFIG
    batch_key: str = cfg.get("batch_key") or smk.config["batch_key"]
    if to_exclude := cfg.get("exclude"):
        adata = adata[~adata.obs["sample"].isin(to_exclude), :]

    old = {}
    for slot in ("obs", "uns", "varp", "obsm"):
        old[slot] = list(getattr(adata, slot).keys())

    key = get_integration_key(smk.params["integration"])
    sc.external.pp.bbknn(adata, batch_key=batch_key)

    if cfg["method"] == "leiden":
        sc.pp.neighbors(adata, use_rep=key)
        sc_utils.sweep_clustering(
            adata,
            lambda x, res, key: sc.tl.leiden(
                x, resolution=res, key_added=key, **cfg["kws"]
            ),
            values=cfg["sweep"],
            prefix="leiden_res",
            distances=adata.obsp["distances"],
        )
    else:
        raise NotImplementedError()
    for slot, key_list in old.items():
        dct = getattr(adata, slot)
        for key in key_list:
            del dct[key]

    del adata.X
    del adata.var
    adata.write_h5ad(smk.output[0])


def prepare_data():
    _ = fn.prepare_data(smk.output[0], smk.output[1], smk.config)


def cellassign():
    adata = ad.read_h5ad(smk.input[0])
    _, markers = fn.get_gs_and_cell_markers(smk.config)
    result = sc_utils.cell_assign_wrapper(
        adata,
        cell_markers=markers,
        model_path=Path(smk.params["model"]),
        **(RCONFIG.get("kws") or {}),
    )
    model = result["model"]
    result["pred"].reset_index().to_csv(smk.output["predictions"], index=False)
    run_metrics = pd.concat([val for val in model.history.values()], axis=1)
    model.save(smk.params["model"], save_anndata=True, overwrite=True)
    run_metrics.reset_index().to_csv(smk.output["run_metrics"], index=False)


def do_dimensionality_reduction():
    cfg = smk.config["DR"]
    adata = ad.read_h5ad(smk.input[0])
    method = smk.params["dr_method"]
    imethod = smk.params["integration_method"]
    kws: dict = cfg["methods"][method].get("kws", {}) or {}
    if imethod == "unintegrated":
        x: np.ndarray = adata.obsm["X_pca"][:, : cfg["n_pcs"]]
        kws["init"] = "normal"
    else:
        x = adata.obsm[get_integration_key(imethod)]
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
