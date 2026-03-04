#!/usr/bin/env ipython
from collections.abc import Callable
from functools import reduce
from pathlib import Path
from typing import Literal

import anndata as ad
import chula_stem.sc_rnaseq as sc_utils
import numpy as np
import pandas as pd
import plotnine as gg
import pymupdf
import scanpy as sc
from chula_stem.plotting import plot_associations
from chula_stem.r_utils import edgeR_wrapper
from pymupdf import Document
from snakemake.script import snakemake as smk

import functions as fn

RNG: int = smk.config["rng"]
RCONFIG: dict = smk.config.get(smk.rule) or {}


# * Utility functions


def get_integration_key(integration_method: str) -> str:
    return f"X_{integration_method}"


def dr_dispatch(
    adata: ad.AnnData, method: str, integration_method: str, kws: dict
) -> Callable[[np.ndarray], np.ndarray]:
    if method.endswith("_tdr"):
        import torchdr as tdr

        if integration_method == "unintegrated":
            x: np.ndarray = adata.obsm["X_pca"][:, : RCONFIG["n_pcs"]]
            kws["init"] = "normal"
        else:
            x = adata.obsm[get_integration_key(integration_method)]
        if method == "umap_tdr":
            dr_obj = tdr.UMAP(**kws)
        elif method == "t-sne_tdr":
            dr_obj = tdr.TSNE(**kws)
        elif method == "pacmap_tdr":
            dr_obj = tdr.PACMAP(**kws)

        else:
            raise ValueError(f"DR method {method} not supported")
        return dr_obj.fit_transform(x)
    else:
        pass
    if integration_method == "unintegrated":
        key = None
    else:
        key = get_integration_key(integration_method)
    if method == "umap":
        neighbor_args = {"n_neighbors", "n_pcs", "knn", "method"}
        keys = list(kws.keys())
        neighbor_kws = {k: kws.pop(k) for k in keys if k in neighbor_args}
        sc.pp.neighbors(adata, use_rep=key, **neighbor_kws)
        sc.tl.umap(adata, **kws)
        return adata.obsm["X_umap"]
    elif method == "t-sne":
        sc.tl.tsne(adata, use_rep=key, use_fast_tsne=True, **kws)
        return adata.obsm["X_tsne"]
    else:
        raise ValueError(f"Method {method} not recognized")


def integrate_data(
    adata: ad.AnnData,
    batch_key,
    name: str,
    extras: dict,
    env,
    params: dict | None = None,
    method: str | None = None,
) -> str:
    """
    Parameters
    ----------
    name : str
        Key in `integration` section of configuration providing parameters for an
        integration method
    extras : dict
        Dictionary for retrieving any intermediate output
    """
    params = params or env["integration"][name]
    method = method or params.get("method", name)
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


def add_all_clusters(
    sample_clusters: list[str],
    adata: ad.AnnData | None = None,
    env: dict | None = None,
    to_category: bool = False,
) -> tuple[pd.DataFrame, ad.AnnData | None]:
    clst_df: pd.DataFrame = reduce(
        lambda x, y: x.merge(y, on="sample", how="outer"),
        [pd.read_csv(f) for f in sample_clusters],
    ).astype("string")
    if to_category:
        clst_df = clst_df.astype("category")
    if adata is not None:
        if env is not None:
            cluster_names = fn.add_clusterings(
                adata, clusters_to_add=env["chosen_clusters"], env=env
            )
        adata.obs = adata.obs.merge(clst_df, on="sample", how="left")
    return clst_df, adata


# * Rules


def gprofiler_enrich():
    sample_level = pd.read_csv(smk.input["sample_level"])
    is_scVI = sample_level["analysis_group"].str.startswith("scVI")
    scVI_proba_threshold: float = RCONFIG.get("scVI_min_proba_de", 0.8)
    scVI_passed_threshold = sample_level["proba_de"] >= scVI_proba_threshold
    sample_level = sample_level[scVI_passed_threshold | ~is_scVI, :]
    sl = fn.profile_de_results(sample_level, level="sample", **RCONFIG)
    sl.to_csv(smk.output["sample_level"])
    tmp = []
    for infile in (Path(p) for p in smk.input["cluster_level"]):
        method = infile.stem.removeprefix("clusters-").removesuffix("_de")
        df = pd.read_csv(infile)
        if "scVI_de" in infile.stem:
            df = df[df["proba_de"] >= scVI_passed_threshold, :]
        cur = fn.profile_de_results(df, level="cluster", **RCONFIG).assign(
            method=method
        )
        tmp.append(cur)
    pd.concat(tmp).to_csv(smk.output["cluster_level"], index=False)


def integrate(adata: ad.AnnData | None = None, cfg: dict | None = None):
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


def cluster_cells(cfg: dict | None = None):
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
    hp_value: int | float = smk.params["hp_value"]
    hp_to_vary = cfg["methods"][method]["vary"][0]
    kws[hp_to_vary] = int(hp_value)
    result: np.ndarray = dr_dispatch(adata, method, integration_method=imethod, kws=kws)
    np.save(smk.output[0], result, allow_pickle=True)


def do_de_samples() -> ad.AnnData | None:
    with open(smk.input["features"], "r") as f:
        feature_idx = f.read().splitlines()
    clst_df: pd.DataFrame = pd.read_csv(smk.input["clusters"]).astype("string")
    adata = ad.read_h5ad(smk.input["adata"])
    adata = adata[:, adata.var.index.isin(feature_idx)]
    adata.obs = adata.obs.merge(clst_df, on="sample", how="left")

    de_spec = smk.params["spec"]
    params = smk.config["do_de_samples"][de_spec]
    method = params.get("method", de_spec)
    query = params.get("query")
    kws = params.get("kws") or {}
    if query:
        kept_ids = adata.obs.query(query).index
        adata = adata[adata.obs.index.isin(kept_ids), :]
    if method == "scVI":
        import scvi

        model: scvi.model.SCVI = scvi.model.SCVI.load(
            smk.input["model"], adata=adata.copy()
        )
        result = (
            model.differential_expression(mode="change", **kws)
            .reset_index(names="gene")
            .rename(columns={"comparison": "contrast"})
        )
    elif method == "edgeR" and (group := kws.pop("group", None)):
        agg_obs = adata.obs.groupby("sample").agg("first")
        adata = sc.get.aggregate(adata, by="sample", func="sum")
        adata.obs = agg_obs
        adata.X = adata.layers["sum"]
        num_de, result = edgeR_wrapper(adata, group=group, **kws)
    else:
        raise NotImplementedError()
    result.to_csv(smk.output[0], index=False)


def train_scvi_permissive(
    adata_file: str | None = None,
    feature_idx_file: str | None = None,
    params=RCONFIG,
    env=smk.config,
):
    """Helper function for training scVI model after reducing adata to the
    permissive feature set
    Input anndata object is `passed_qc`
    """
    import scvi

    adata_file = adata_file or smk.input[0]
    feature_idx_file = feature_idx_file or smk.input[1]

    adata = ad.read_h5ad(adata_file)
    with open(feature_idx_file, "r") as f:
        features = f.read().splitlines()
    adata = adata[:, adata.var.index.isin(features)].copy()
    batch_key = env["batch_key"]
    tmp = {}
    integrate_data(adata, batch_key, "", tmp, env, params=params, method="scVI")
    model: scvi.model.SCVI = tmp["model"]
    model.save(smk.output["path"], save_anndata=False, overwrite=True)


def save_other_dotplots():
    adata = ad.read_h5ad(smk.input["adata"])
    cols_plot = ["cfs_community", "cellassign_prediction", "patient", "sample"]
    fn.add_cfs_community(
        adata,
        old_clusters=Path(smk.input["ind_clustering"]),
        cfs_community_file=Path(smk.input["communities"]),
        column="cfs_community",
    )
    cellassign_predictions: pd.DataFrame = pd.read_csv(smk.input["predictions"])
    cellassign_metrics: pd.DataFrame = (
        pd.read_csv(smk.input["run_metrics"])
        .melt(id_vars="epoch")
        .rename(columns={"variable": "metric"})
    )
    metrics_plot = (
        gg.ggplot(cellassign_metrics, gg.aes(x="epoch", y="value"))
        + gg.geom_line()
        + gg.facet_wrap("metric", scales="free", ncol=2)
        + gg.ggtitle("Cellassign training metrics")
        + gg.theme(figure_size=(15, 20))
    )
    metrics_plot.save(smk.output[1])

    adata.obs = adata.obs.merge(
        cellassign_predictions.loc[:, ["index", "PREDICTION"]],
        left_index=True,
        right_on="index",
        how="left",
    ).rename(columns={"PREDICTION": "cellassign_prediction"})
    fn.make_cluster_dotplots(
        adata,
        filename=smk.output[0],
        markers={},
        env=smk.config,
        additional_groups=cols_plot,
        group_rotation=90,
    )


def enrich_clusters():
    adata = ad.read_h5ad(smk.input[0])
    fn.enrich_clusters(
        adata,
        cfg=RCONFIG,
        env=smk.config,
        gs_out=smk.output[0],
        marker_out=smk.output[1],
    )


def do_de_clusters():
    method = smk.params["method"]
    cfg = RCONFIG.get(method) or {}
    adata = ad.read_h5ad(smk.input["adata"])
    res = fn.do_de_clusters(
        method,
        adata,
        cfg,
        smk.config,
        model_file=smk.input["model"],
        features=smk.input["features"],
    )
    if method == "edgeR":
        res["de_counts"].to_csv(
            f"{smk.params['outdir']}/edgeR_de_gene_counts.csv", index=False
        )
    res["top_de"].to_csv(smk.output[0], index=False)


def save_sample_dotplots():
    adata = ad.read_h5ad(smk.input["adata"])
    clst_df, adata = add_all_clusters(smk.input["clusters"], adata)
    assert isinstance(adata, ad.AnnData)
    fn.make_cluster_dotplots(
        adata,
        filename=smk.output[0],
        markers={},
        env=smk.config,
        additional_groups=[c for c in clst_df.columns if c != "sample"],
        group_rotation=90,
        with_samples=False,
    )


def gather_sample_clusters():
    clst_df, _ = add_all_clusters(smk.input["clusters"])
    clst_df.to_csv(smk.output[0], index=False)
    doc: Document = pymupdf.open()
    for d in smk.input["plots"]:
        doc.insert_file(d)
    doc.save(smk.output[2])


def find_gene_associations():
    adata = ad.read_h5ad(smk.input["adata"])
    features: list[str] = Path(smk.input["features"]).read_text().splitlines()
    _, adata = add_all_clusters(
        smk.input["clusters"], adata, smk.config, to_category=True
    )
    assert isinstance(adata, ad.AnnData)
    outdir: Path = Path(smk.params["outdir"])
    kws: dict = RCONFIG["kws"]
    spec: dict | None = None
    name = smk.params["name"]
    for spec_tmp in RCONFIG["spec"]:
        if name == spec_tmp["name"]:
            spec = spec_tmp
    if spec is None:
        raise ValueError("Spec is none for some reason")
    outfile = smk.output[1]
    qgenes: list = spec["genes"]
    if query := spec.get("query"):
        kept: pd.DataFrame = adata.obs.query(query)
        filtered: ad.AnnData = adata[adata.obs_names.isin(kept.index), :]
    else:
        filtered = adata
    assert filtered.shape[0] > 0
    if Path(smk.output[0]).exists():
        result = pd.read_csv(smk.output[0])
    else:
        result: pd.DataFrame = sc_utils.find_proportional_genes_rho(
            filtered, qgenes, **kws
        )
        result.to_csv(smk.output[0], index=False)
    plots = plot_associations(
        adata,
        groupby=spec.get("group_by", "sample"),
        assoc_df=result,
        n=RCONFIG.get("keep", 10),
        style=RCONFIG.get("style", "heatmap"),
    )
    doc: Document = pymupdf.open()
    for i, plot in enumerate(plots):
        tmp = outdir / f"{i}.pdf"
        plot.savefig(tmp, bbox_inches="tight")
        doc.insert_file(tmp)
        tmp.unlink()
    doc.save(outfile)


# * Entry
if rule_fn := globals().get(smk.rule):
    rule_fn()
elif smk.rule.startswith("do_dimensionality_reduction"):
    do_dimensionality_reduction()
