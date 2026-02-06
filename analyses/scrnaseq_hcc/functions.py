#!/usr/bin/env python

from pathlib import Path

import anndata as ad
import polars as pl
import scanpy as sc
from chula_stem.sc_rnaseq import annotate_adata_vars, distance_by_mads
from chula_stem.utils import read_existing
from loguru import logger


def prepare_data(env: dict) -> ad.AnnData:
    def helper(file):
        outdir = Path(env["outdir"])
        tmp = []
        manifest = pl.read_csv(env["manifest"])
        target = "filtered_feature_bc_matrix.h5"
        for row in manifest.iter_rows(named=True):
            patient = row["sample_name"]
            stype = row["type"]
            suffix = f"{row['type']}-{row['treatment']}"
            data_path: Path = (
                outdir / patient / "processed" / f"cellranger_{suffix}/{target}"
            )
            if not data_path.exists():
                logger.warning(
                    "The files for sample {} ({}) don't exist", patient, suffix
                )
                continue
            current: ad.AnnData = sc.read_10x_h5(data_path)
            current.obs_names_make_unique()
            current.var_names_make_unique()
            current.obs.loc[:, "sample"] = f"{patient}_{suffix}"
            current.obs.loc[:, "patient"] = patient
            current.obs.loc[:, "type"] = stype
            current.obs.loc[:, "treatment"] = row["treatment"]
            tmp.append(current)
        if extra := env["files"].get("extras"):
            for key, path in extra.items():
                patient, stype = key.split("|")
                if not Path(path).exists():
                    logger.warning("Path to extra file {} does not exist", path)
                    continue
                current = ad.read_h5ad(path)
                current.obs.loc[:, "sample"] = patient
                current.obs.loc[:, "type"] = stype
                tmp.append(current)
        adata = ad.concat(tmp, merge="first", index_unique="-")
        annotate_adata_vars(adata, "gene_ids", savepath=Path(env["gene_reference"]))
        print(adata.var.columns)
        sc.pp.calculate_qc_metrics(adata, inplace=True, qc_vars=["mito"])
        # TODO: should you change the grouping col?
        distance_by_mads(
            adata,
            ["total_counts", "n_genes_by_counts", "pct_counts_mito"],
            group_keys=None,
            inplace=True,
        )
        print(adata.uns)
        adata.write_h5ad(file)
        adata = ad.read_h5ad(file, backed=True)
        return adata

    adata = read_existing(
        Path(env["files"]["combined"]), helper, lambda x: ad.read_h5ad(x, backed=True)
    )
    print(f"Samples present: {set(adata.obs["sample"])}")
    return adata


# TODO: write a function to add mads calculation and mitochondrial gene annotation
