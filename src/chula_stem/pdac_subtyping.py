#!/usr/bin/env ipython
from collections.abc import Sequence
from pathlib import Path
from typing import Callable

import anndata as ad
import click
import numpy as np
import pandas as pd
import pandera.pandas as pa
from awkward import count
from beartype import beartype
from jax.numpy import single
from scipy import sparse
from scipy.special import expit
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import train_test_split

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
    adata.X = sparse.csc_array(adata.X)
    return adata


def group_by_index(files: Sequence[str], index_col: str, **kws) -> list[list]:
    tmp = {"index": [], "file": []}
    for f in files:
        cur = pd.read_csv(f, **kws)
        tmp["index"].append(tuple(cur[index_col]))
        tmp["file"].append(f)
    return [g["file"].tolist() for _, g in pd.DataFrame(tmp).groupby("index")]


class MoffittBasal:
    def __init__(self, threshold: float = 0.5) -> None:
        self.threshold: float = threshold
        self.intercept = -7.16
        self.adata: ad.AnnData | None = None
        self.genes: pd.Index | None = None
        self.x: np.ndarray | None = None
        self.df: pd.DataFrame | None = None

    def fit(
        self, adata: ad.AnnData, gene_col: str | None = None, layer: str | None = None
    ):
        # Was CTSL2 in the original paper
        self.adata = adata
        ctsl = "CTSL" if "CTSL" in adata.var_names else "CTSL2"
        self.genes = pd.Index(adata.var[gene_col]) if gene_col else adata.var_names
        data = {
            "A": [
                "CD109",
                "SLC2A1",
                "KRT16",
                ctsl,
                "KRT6A",
                "B3GNT5",
                "MET",
                "CHST6",
                "SERPINB5",
                "DCBLD2",
                "IL20RB",
                "PPP1R14C",
                "NAB1",
                "MSLN",
            ],
            "B": [
                "GPR160",
                "AGR2",
                "SLC44A4",
                "TMEM45B",
                "BCAS1",
                "VSIG2",
                "TFF3",
                "PLA2G10",
                "HPGD",
                "PLS1",
                "FAM3D",
                "SYTL2",
                "PLEKHA6",
                "CAPN9",
            ],
            "coeff": [
                0.87,
                1.22,
                0.52,
                1.43,
                0.70,
                0.41,
                0.72,
                0.80,
                0.76,
                1.40,
                1.33,
                1.58,
                0.41,
                1.58,
            ],
        }
        self.df = pd.DataFrame(data)
        for col in ("A", "B"):
            if not self.df[col].isin(self.genes).all():
                missing = self.df[col][~self.df[col].isin(self.genes)]
                raise ValueError(
                    f"The following {col} genes are missing from adata: {missing.tolist()}"
                )
        # "Estimated Increase in odds of basal class membership when A > B"
        self.df.loc[:, "odds_increase"] = np.exp(self.df["coeff"])
        a_locs = [self.genes.get_loc(g) for g in self.df["A"]]
        b_locs = [self.genes.get_loc(g) for g in self.df["B"]]
        assert isinstance(self.adata.var, pd.DataFrame)
        expr: np.ndarray = self.adata.X if layer is None else self.adata.layers[layer]

        if not isinstance(expr, np.ndarray) and "toarray" in dir(expr):
            expr = expr.toarray()
        a_vals = expr[:, a_locs]
        b_vals = expr[:, b_locs]
        self.x = a_vals > b_vals

    def predict(
        self, inplace: bool = True, _x: np.ndarray | None = None
    ) -> np.ndarray | None:
        if self.x is None or self.df is None:
            raise ValueError("Must call `fit` first")
        x = self.x if _x is None else _x
        power = np.matmul(x, self.df["coeff"].values.reshape((-1, 1)))
        proba = expit(self.intercept + power)
        odds = proba / (1 - proba)
        # odds2 = np.log(np.matmul(mask, df["odds_increase"].values.reshape((-1, 1))))
        if not inplace:
            return proba
        self.adata.obs.loc[:, "moffitt_pr_basal"] = proba
        self.adata.obs.loc[:, "moffitt_is_basal"] = np.where(
            proba > self.threshold, [1], [0]
        )

    def tune(
        self,
        y_true: str,
        metric_fn: Callable = balanced_accuracy_score,
        n: int = 5,
        **kws,
    ) -> None:
        """Use sklearn's TunedThresholdClassifier to choose the best cut-off for class labels"""
        if self.x is None or self.df is None:
            raise ValueError("Must call `fit` first")
        y_true = self.adata.obs[y_true].values

        bests = []
        for i in range(n):
            x_train, x_test, y_train, y_test = train_test_split(
                self.x, y_true, stratify=y_true, **kws
            )
            n_train = len(x_train)

            def call(t: float):
                proba = self.predict(_x=x_train, inplace=False)
                pred = np.where(proba > t, [1], [0])
                return metric_fn(y_train, pred)

            thresholds = np.arange(0, 1, 0.05)
            scores = np.array([call(t) for t in thresholds])
            print("Scores at thresholds")
            print(
                pd.DataFrame({"Threshold": thresholds, "Scores": np.round(scores, 3)})
            )
            best = thresholds[np.argmax(scores)]

            x_pred = np.where(
                self.predict(_x=x_test, inplace=False) > best,
                [1],
                [0],
            )
            print(f"Score on test set: {metric_fn(y_test, x_pred)}")
            bests.append(best)

        self.threshold = np.mean(bests)


            "CD109",
            "SLC2A1",
            "KRT16",
            "CTSL",  # Was originally CTSL2
            "KRT6A",
            "B3GNT5",
            "MET",
            "CHST6",
            "SERPINB5",
            "DCBLD2",
            "IL20RB",
            "PPP1R14C",
            "NAB1",
            "MSLN",
        ],
        "B": [
            "GPR160",
            "AGR2",
            "SLC44A4",
            "TMEM45B",
            "BCAS1",
            "VSIG2",
            "TFF3",
            "PLA2G10",
            "HPGD",
            "PLS1",
            "FAM3D",
            "SYTL2",
            "PLEKHA6",
            "CAPN9",
        ],
        "coeff": [
            0.87,
            1.22,
            0.52,
            1.43,
            0.70,
            0.41,
            0.72,
            0.80,
            0.76,
            1.40,
            1.33,
            1.58,
            0.41,
            1.58,
        ],
    }
    assert isinstance(adata.var, pd.DataFrame)
    genes: pd.Series = adata.var[gene_col] if gene_col else pd.Series(adata.var_names)
    x: np.ndarray = adata.X if layer is None else adata.layers[layer]
    if not isinstance(x, np.ndarray) and "toarray" in dir(x):
        x = x.toarray()

    df = pd.DataFrame(data)
    for col in ("A", "B"):
        if not df[col].isin(genes).all():
            missing = df[col][~df[col].isin(genes)]
            raise ValueError(
                f"The following {col} genes are missing from adata: {missing.tolist()}"
            )
    # "Estimated Increase in odds of basal class membership when A > B"
    df.loc[:, "odds_increase"] = np.exp(df["coeff"])

    a_vals = x[:, genes.isin(df["A"])]
    b_vals = x[:, genes.isin(df["B"])]
    mask = a_vals > b_vals
    power = np.matmul(mask, df["odds_increase"].values.reshape((-1, 1)))
    proba = expit(intercept + power)
    # print(odds_basal)
    # proba = (odds_basal / (1 + odds_basal)).flatten()
    if not inplace:
        return proba
    adata.obs.loc[:, "moffitt_pr_basal"] = proba
