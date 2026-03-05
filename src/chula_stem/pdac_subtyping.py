#!/usr/bin/env ipython
from collections.abc import Sequence
from pathlib import Path

import anndata as ad
import click
import numpy as np
import pandas as pd
import pandera.pandas as pa
from awkward import count
from beartype import beartype

from chula_stem.r_utils import tximport

"[2026-03-04 Wed] Interested in GATA6 and RAS want to score a transcriptome based on activity of these genes"

# TODO: [2026-03-04 Wed] for now also just plot GATA6 expression distribution
# among samples just do heatmap

from chula_stem.utils import get_ensembl_gene_data

SCHEMA: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "file": pa.Column(
            str, pa.Check(lambda x: Path(x).exists(), element_wise=True), unique=True
        ),
        "sample": pa.Column(str, nullable=False, unique=False),
        "type": pa.Column(
            str, pa.Check(lambda x: x.isin({"SV", "salmon", "kallisto"}))
        ),
    },
    unique=["sample", "type"],
)


@beartype
def read_manifest(manifest: pd.DataFrame, tx2gene: pd.DataFrame) -> ad.AnnData:
    SCHEMA.validate(manifest)

    sv_check: set = {"SV"}
    no_expr: list[str] = (
        manifest.groupby("sample")
        .agg({"type": set})
        .query("type == @sv_check")
        .index.to_list()
    )
    count_adatas: list[ad.AnnData] = []
    for group, df in manifest.groupby("type"):
        names = df["sample"]
        files = df["file"]
        if group in {"salmon", "kallisto"}:
            index_col = "Name" if group == "salmon" else "target_id"
            index_groups = group_by_index(files, index_col, sep="\t")
            # necessary to avoid errors with files not having the same index
            for g in index_groups:
                cur = df.query("file.isin(@g)")
                adata = tximport(
                    files=list(cur["file"]),
                    tx2gene=tx2gene,
                    sample_names=list(cur["sample"]),
                    type=group,
                )
                adata.obs.loc[:, "rnaseq_available"] = ~adata.obs.index.isin(no_expr)
                count_adatas.append(adata)
        elif group == "SV":
            # TODO: want to read in the SV data as a dataframe and add it to adata.obs
            pass
    adata = ad.concat(count_adatas, merge="first", join="outer", uns_merge="same")
    if no_expr:
        dummy = ad.AnnData(
            X=np.zeros((len(no_expr), adata.shape[1])),
            var=pd.DataFrame(index=adata.var_names),
            obs=pd.DataFrame({"rnaseq_available": False}, index=no_expr),
        )
        adata = ad.concat([adata, dummy], merge="first", join="outer", uns_merge="same")
    return adata


def group_by_index(files: Sequence[str], index_col: str, **kws) -> list[list]:
    tmp = {"index": [], "file": []}
    for f in files:
        cur = pd.read_csv(f, **kws)
        tmp["index"].append(tuple(cur[index_col]))
        tmp["file"].append(f)
    return [g["file"].tolist() for _, g in pd.DataFrame(tmp).groupby("index")]
