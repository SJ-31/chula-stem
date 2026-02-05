#!/usr/bin/env ipython

from pathlib import Path

import anndata as ad
import polars as pl
import scanpy as sc
from chula_stem.utils import read_existing
from loguru import logger


def prepare_data(env: dict) -> ad.AnnData:
    def helper(file):
        outdir = Path(env["outdir"])
        tmp = []
        manifest = pl.read_csv(env["manifest"])
        target = "filtered_feature_bc_matrix.h5"
        for row in manifest.iter_rows(named=True):
            sample = row["sample_name"]
            stype = row["type"]
            data_path: Path = (
                outdir / sample / "processed" / f"cellranger_{stype}/{target}"
            )
            if not data_path.exists():
                logger.warning(
                    "The files for sample {} ({}) don't exist", sample, stype
                )
                continue
            current: ad.AnnData = sc.read_10x_h5(data_path)
            current.obs_names_make_unique()
            current.var_names_make_unique()
            current.obs.loc[:, "sample"] = sample
            current.obs.loc[:, "type"] = stype
            tmp.append(current)
        if extra := env["files"].get("extras"):
            for key, path in extra.items():
                sample, stype = key.split("|")
                if not Path(path).exists():
                    logger.warning("Path to extra file {} does not exist", path)
                    continue
                current = ad.read_h5ad(path)
                current.obs.loc[:, "sample"] = sample
                current.obs.loc[:, "type"] = stype
                tmp.append(current)
        adata = ad.concat(tmp)
        adata.write_h5ad(file)
        return adata

    adata = read_existing(Path(env["files"]["combined"]), helper, ad.read_h5ad)
    print(f"Samples present: {set(adata.obs["sample"])}")
    return adata
