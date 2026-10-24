#!/usr/bin/env ipython

import faulthandler
from collections.abc import Callable
from pathlib import Path

import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd
import yte
from chula_stem.r_utils import tmm
from loguru import logger
from pyhere import here
from sklearn.metrics import adjusted_rand_score

with open(here("analyses", "cluster_reproducibility", "env.yaml")) as f:
    env = yte.process_yaml(f)

faulthandler.enable()


def decoupler_wrapper(
    fn: Callable, **kws
) -> Callable[[ad.AnnData, pd.DataFrame], np.ndarray]:
    def inner(adata: ad.AnnData, net: pd.DataFrame) -> np.ndarray:
        logger.debug("Got into inner")
        adata.obsm = {}
        logger.debug("Starting fn...")
        fn(adata, net, **kws)
        key_added = [k for k in adata.obsm.keys() if k.startswith("score_")][0]
        acts: pd.DataFrame = adata.obsm[key_added]
        logger.info("Activity values: {}", acts)
        return acts.idxmax(axis=1).values

    return inner


def score_cluster_reproduction(
    adata: ad.AnnData,
    gs_scorer: Callable[[ad.AnnData, pd.DataFrame], np.ndarray],
    cluster_col: str = "best_cc",
    n_per_set: int = 30,
):
    """Evaluate how well gene-set scoring function `gs_scorer` is able to reproduce
    the clustering that gave rise to its gene sets

    Parameters
    ----------
    n_per_set : int
        Number of DE genes to use when characterizing each cluster

    """
    de_df: pd.DataFrame = adata.uns["top_de"]
    de_df = de_df.assign(
        contrast=de_df["contrast"]
        .str.replace(" vs.*", "", regex=True)
        .str.removeprefix(cluster_col)
    )
    if de_df.empty:
        logger.warning("No DE genes available between clusters")
        return
    gene_sets = {
        g: df.sort_values("logFC", ascending=False).head(n_per_set)["gene"]
        for g, df in de_df.groupby("contrast")
    }
    net = pd.DataFrame(
        {"source": gene_sets.keys(), "target": gene_sets.values()}
    ).explode("target")
    replication: np.ndarray = gs_scorer(adata, net)
    true = adata.obs[cluster_col]
    return adjusted_rand_score(labels_true=true, labels_pred=replication)


workdir = here(*env["workdir"])
files: Path = here(*env["outs"]["clustered"])

methods = {
    "gsva": decoupler_wrapper(dc.mt.gsva, tmin=3),
    "ora": decoupler_wrapper(dc.mt.ora, tmin=3),
    "zscore": decoupler_wrapper(dc.mt.zscore, tmin=3),
    "aucell": decoupler_wrapper(dc.mt.aucell, tmin=3),
    "ulm": decoupler_wrapper(dc.mt.ulm, tmin=3),
}

results = {"dataset": [], "method": [], "ari": []}
for file in files.glob("*.h5ad"):
    adata = ad.read_h5ad(file)
    adata.X = tmm(adata, log=True)
    logger.info("Current adata {}", adata)
    logger.info("DE genes {}", adata.uns["top_de"].iloc[:, :3])
    for k, v in methods.items():
        results["dataset"].append(file.stem)
        results["method"].append(k)
        logger.info("Scoring {} with {}", file.stem, k)
        results["ari"].append(score_cluster_reproduction(adata, gs_scorer=v))
df = pd.DataFrame(results)
df.to_csv(workdir / "result.csv", index=False)
