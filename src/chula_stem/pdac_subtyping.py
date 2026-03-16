#!/usr/bin/env ipython
from collections.abc import Sequence
from pathlib import Path
from typing import Callable, Literal

import anndata as ad
import click
import numpy as np
import pandas as pd
import pandera.pandas as pa
import pandera.polars as pal
import polars as pl
from beartype import beartype
from scipy import sparse
from scipy.special import expit
from scipy.stats import expon, goodness_of_fit, iqr
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import train_test_split

from chula_stem.r_utils import tximport

"[2026-03-04 Wed] Interested in GATA6 and RAS want to score a transcriptome based on activity of these genes"

# TODO: [2026-03-04 Wed] for now also just plot GATA6 expression distribution
# among samples just do heatmap


SCHEMA: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "file": pa.Column(
            str, pa.Check(lambda x: Path(x).exists(), element_wise=True), unique=True
        ),
        "sample": pa.Column(str, nullable=False, unique=False),
        "type": pa.Column(
            str,
            pa.Check(lambda x: x.isin({"SV", "salmon", "kallisto", "purecn"})),
        ),
    },
    unique=["sample", "type"],
)


# TODO: also wanna collect the PureCN summary files and plot the comments (as a heatmap),
# and purity + ploidy of the samples


@beartype
def read_manifest(
    manifest: pd.DataFrame, tx2gene: pd.DataFrame, purecn_prefix: str = ""
) -> ad.AnnData:
    SCHEMA.validate(manifest)

    no_expr: list[str] = (
        pl.from_pandas(manifest)
        .group_by("sample")
        .agg(pl.col("type"))
        .filter(
            pl.col("type").list.set_intersection(["kallisto", "salmon"]).list.len() == 0
        )["sample"]
        .to_list()
    )

    adata: ad.AnnData = get_counts(manifest, tx2gene, no_expr)
    if no_expr:
        dummy = ad.AnnData(
            X=np.zeros((len(no_expr), adata.shape[1])),
            var=pd.DataFrame(index=adata.var_names),
            obs=pd.DataFrame({"rnaseq_available": False}, index=no_expr),
        )
        adata = ad.concat([adata, dummy], merge="first", join="outer", uns_merge="same")

    adata.X = sparse.csc_array(adata.X)
    return adata


def add_kras_imbalance(
    manifest: pd.DataFrame, adata: ad.AnnData, purecn_prefix: str = ""
) -> None:
    purecn_df = manifest.query("type == 'purecn'")
    if not purecn_df.empty:
        if "KRAS_imbalance_subtype" in adata.obs:
            adata.obs = adata.obs.drop("KRAS_imbalance_subtype", axis="columns")
        kras_scores = call_kras_imbalance(purecn_df, purecn_prefix=purecn_prefix)
        adata.obs = adata.obs.merge(kras_scores["subtype"], on="sample", how="left")
        for key in ("kras_variants", "kras_genes", "summary"):
            df = kras_scores[key]
            if not df.empty:
                adata.uns[f"purecn_{key}"] = df
            else:
                adata.uns[f"purecn_{key}"] = pd.DataFrame()
    else:
        print("Warning: No samples with PureCN directories in manifest")


def group_by_index(files: Sequence[str], index_col: str, **kws) -> list[list]:
    tmp = {"index": [], "file": []}
    for f in files:
        cur = pd.read_csv(f, **kws)
        tmp["index"].append(tuple(cur[index_col]))
        tmp["file"].append(f)
    return [g["file"].tolist() for _, g in pd.DataFrame(tmp).groupby("index")]


def get_counts(
    manifest: pd.DataFrame, tx2gene: pd.DataFrame, no_expr: list[str]
) -> ad.AnnData:
    count_adatas: list[ad.AnnData] = []
    for group, df in manifest.groupby("type"):
        files = df["file"]
        if group in {"salmon", "kallisto"}:
            index_col = "Name" if group == "salmon" else "target_id"
            index_groups = group_by_index(files, index_col, sep="\t")
            # necessary to avoid errors with files not having the same index
            for index_group in index_groups:
                cur = df.loc[df["file"].isin(index_group), :]
                adata = tximport(
                    files=list(cur["file"]),
                    tx2gene=tx2gene,
                    sample_names=list(cur["sample"]),
                    type=group,
                )
                adata.obs.loc[:, "rnaseq_available"] = ~adata.obs.index.isin(no_expr)
                count_adatas.append(adata)
    adata = ad.concat(count_adatas, merge="first", join="outer", uns_merge="same")
    return adata


# * Moffitt


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
        self.adata.obs.loc[:, "moffitt_is_basal"] = self.adata.obs[
            "moffitt_is_basal"
        ].astype("string")

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


# * KRAS instability


def instability_score_one(
    v_df: pd.DataFrame,
    g_df: pd.DataFrame,
) -> Literal["wt", "balanced", "minor", "major", "NA"]:
    """
    Provide a KRAS imbalance score for the given sample, following
    https://www.nature.com/articles/s41588-019-0566-9/figures/13

    The `major` allele is in the notes and variable names refers to the wild-type
    reference allele (as the calls were filtered to be somatic only)

    Notes
    -----
    The samples fail classification (assigned 'NA') in any of the following situations
        The KRAS gene was deleted completely (total CN == 0)
        TODO: find other situations of this
    """
    v_df = v_df.assign(major_c=v_df["ML.C"] - v_df["ML.M.SEGMENT"])
    g_df = g_df.assign(majo_C=g_df["C"] - g_df["M"])
    if v_df.empty and g_df.empty:
        return "wt"
    elif v_df.empty:
        if g_df.shape[0] != 1:
            print(g_df)
            raise ValueError("Shape of gene_df after filtering should be 1")
        total_cn = g_df["C"].item()
        minor_cn = g_df["M"].item()
        major_cn = total_cn - minor_cn
        # If loh, the gene bearing variants was lost
        if total_cn == 0:  # Unsure how to categorize deletions
            return "NA"
        if major_cn == minor_cn:
            return "balanced"
        if major_cn > minor_cn:  # TODO: is assuming that loh is wt correct?
            return "wt"
        return "NA"

    n_vars = v_df.shape[0]
    m_cn: pd.Series = v_df["ML.M.SEGMENT"]
    balance_count = m_cn == v_df["major_C"]
    if (balance_count.sum() / n_vars) > 0.5:
        return "balanced"

    no_wt: pd.Series = v_df["major_C"] == 0
    wt_one: pd.Series = v_df["major_C"] == 1

    minor_ins_count = (no_wt & m_cn == 1) | (wt_one & m_cn == 2)
    if (minor_ins_count.sum() / n_vars) > 0.5:
        return "minor"
    major_ins_count = (no_wt & (m_cn > 1)) | (wt_one & (m_cn > 2))
    if (major_ins_count.sum() / n_vars) > 0.5:
        return "major"
    wt_greater = v_df["major_C"] > m_cn
    if (wt_greater.sum() / n_vars) > 0.5:
        return "wt"  # TODO: what about situation when WT copies > mutant copies?
    return "NA"


GENE_DF_DTYPES: dict = {
    "chr": "string",
    "type": "string",
    "loh": "boolean",
    "M.flagged": "boolean",
    "C.flagged": "boolean",
    "M": "Int16",
}
VARIANT_DF_DTYPES: dict = {
    "gene.symbol": "string",
    "POSTERIOR.SOMATIC": "Float64",
    "M.SEGMENT.FLAGGED": "boolean",
    "ML.M.SEGMENT": "Int64",
    "ML.C": "Float64",
}


@beartype
def call_kras_imbalance(
    manifest: pd.DataFrame,
    purecn_prefix: str = "",
) -> dict[str, pd.DataFrame]:
    results: dict = {}
    tmp_subtype = {"sample": [], "KRAS_imbalance_subtype": []}
    gene_schema = pa.DataFrameSchema(
        {
            "gene.symbol": pa.Column(str),
            "loh": pa.Column("boolean", nullable=True),
            "M": pa.Column("Int16", nullable=True),
            "M.flagged": pa.Column("boolean", nullable=True),
            "C.flagged": pa.Column("boolean", nullable=True),
        }
    )
    variant_schema = pa.DataFrameSchema(
        {
            "gene.symbol": pa.Column("string"),
            "POSTERIOR.SOMATIC": pa.Column(
                "Float64", nullable=True
            ),  # To check whether variant is somatic
            "M.SEGMENT.FLAGGED": pa.Column(
                "boolean", nullable=True
            ),  # To ignore unreliable segments
            "ML.M.SEGMENT": pa.Column("Int64", nullable=True),  # Minor copy number
            "ML.C": pa.Column("Float64", nullable=True),  # Total copy number
        },
        drop_invalid_rows=True,
    )

    df = manifest.pivot(columns="type", index="sample", values="file")
    variant_dfs, gene_dfs, summary_dfs = [], [], []
    for index, series in df.iterrows():
        assert isinstance(index, str)
        dir = Path(series["purecn"])  # output of predictSomatic
        assert dir.exists() and dir.is_dir()
        v_result = dir / f"{purecn_prefix}{index}_variants.csv"
        g_result = (
            dir / f"{purecn_prefix}{index}_genes.csv"
        )  # output of callAlterations
        summary = dir / f"{purecn_prefix}{index}.csv"
        for file in (v_result, g_result, summary):
            if not file.exists():
                raise ValueError(f"PureCN file {file} is missing")
        summary_dfs.append(pd.read_csv(summary, dtype={"Comment": "string"}))

        tmp_subtype["sample"].append(index)
        v_df: pd.DataFrame = pd.read_csv(v_result, dtype=VARIANT_DF_DTYPES)
        g_df: pd.DataFrame = pd.read_csv(g_result, dtype=GENE_DF_DTYPES)

        gene_schema.validate(g_df)
        variant_schema.validate(v_df, lazy=True)
        sym_query = "`gene.symbol` == 'KRAS'"
        v_df = v_df.query(
            f"{sym_query} & (`POSTERIOR.SOMATIC` >= 0.8) & ~`M.SEGMENT.FLAGGED`"
        )
        # Cutoff for somatic calls provided by PureCN manual
        g_df = g_df.query(f"{sym_query} & ~`M.flagged` & ~`C.flagged`")
        variant_dfs.append(v_df.assign(sample=index))
        gene_dfs.append(g_df.assign(sample=index))

        tmp_subtype["KRAS_imbalance_subtype"].append(
            instability_score_one(v_df=v_df, g_df=g_df)
        )
    results["kras_genes"] = pd.concat(gene_dfs)
    results["summary"] = pd.concat(summary_dfs)
    results["subtype"] = pd.DataFrame(tmp_subtype)
    results["kras_variants"] = pd.concat(variant_dfs)
    return results


# * Waddell

# TODO: want to visualize these data as well

CHR_LENGTHS: dict = {
    "1": 248956422,
    "2": 242193529,
    "3": 198295559,
    "4": 190214555,
    "5": 181538259,
    "6": 170805979,
    "7": 159345973,
    "8": 145138636,
    "9": 138394717,
    "10": 133797422,
    "11": 135086622,
    "12": 133275309,
    "13": 114364328,
    "14": 107043718,
    "15": 101991189,
    "16": 90338345,
    "17": 83257441,
    "18": 80373285,
    "19": 58617616,
    "20": 64444167,
    "21": 46709983,
    "22": 50818468,
    "X": 156040895,
    "Y": 57227415,
}


def add_waddell_classification(manifest: pd.DataFrame, adata: ad.AnnData) -> None:
    sv_schema = pal.DataFrameSchema(
        {
            "Loc": pal.Column("string"),
            "Alt": pal.Column("string"),
            "Ref": pal.Column("string"),
            "SVTYPE": pal.Column("string"),
        }
    )
    has_sv = manifest.query("type == 'SV'")
    if not has_sv.empty:
        if "waddell_subtype" in adata.obs:
            adata.obs = adata.obs.drop("waddell_subtype", axis="columns")
        tmp = {"sample": [], "waddell_subtype": []}
        for _, row in has_sv.iterrows():
            tmp["sample"].append(row["sample"])
            tmp["waddell_subtype"].append(classify_waddell_one(row["file"], sv_schema))
        df = pd.DataFrame(tmp)
        adata.obs = adata.obs.merge(df, on="sample", how="left")
    else:
        print("Warning: No samples with SV data in manifest")


def classify_waddell_one(
    file: str, schema: pal.DataFrameSchema
) -> Literal["stable", "scattered", "unstable", "locally_rearranged", "NA"]:
    df = pl.read_csv(file, separator="\t" if file.endswith(".tsv") else ",")
    if df.is_empty():
        return "stable"
    schema.validate(df)
    formatted = (
        df.unique(["Loc", "Alt", "Ref"])
        .filter(pl.col("SVTYPE") != "INS")
        .with_columns(
            pl.col("Loc").str.split_exact(":", 1).struct.rename_fields(["chr", "pos"])
        )
        .unnest("Loc")
        .with_columns(pl.col("pos").cast(pl.Int64))
    )
    if formatted["chr"].str.starts_with("chr").all():
        formatted = formatted.with_columns(pl.col("chr").str.strip_prefix("chr"))
    sv_total: int = formatted.height
    c_df: pl.DataFrame = check_sv_clustering_significance(formatted)
    l_df: pl.DataFrame = check_locally_rearranged(formatted)
    if sv_total < 50 and not c_df["clustering_significant"].any():
        return "stable"
    if 50 <= sv_total < 200 and not c_df["clustering_significant"].any():
        return "scattered"
    if sv_total >= 200 and not c_df["clustering_significant"].any():
        return "unstable"
    if sv_total > 50 and l_df["locally_rearranged"].any():
        return "locally_rearranged"
    return "NA"


def check_locally_rearranged(sv_df: pl.DataFrame) -> pl.DataFrame:
    """
    Identify locally rearranged chromosomes, according to [1]

    1. Calculate SV rate per MB for each chromosome
    2. Mark chromosomes are rearranged if their SV rate > 5 *IQR(SV rate) + quantile(SV rate)

    References
    ----------
    [1] Waddell, N., Pajic, M., Patch, AM. et al. Whole genomes redefine the mutational landscape of pancreatic cancer. Nature 518, 495–501 (2015). https://doi.org/10.1038/nature14169
    """
    grouped = (
        sv_df.group_by("chr")
        .agg(pl.col("pos").len().alias("n"))
        .with_columns(pl.col("chr").replace_strict(CHR_LENGTHS).alias("length"))
        .with_columns((pl.col("n") / pl.col("length") * 10**6).alias("rate_mb"))
    )
    print(grouped["rate_mb"])
    iqr_val = iqr(grouped["rate_mb"])
    upper_quartile = np.quantile(grouped["rate_mb"], 0.75)
    grouped = grouped.with_columns(
        (pl.col("rate_mb") > 5 * iqr_val + upper_quartile).alias("locally_rearranged")
    )
    return grouped


def check_sv_clustering_significance(
    sv_df: pl.DataFrame, threshold: float = 0.0001
) -> pl.DataFrame:
    """
    Method adapted from Box 2, "Clustering of breakpoints" of [1]

    For each chromosome,

    1. Order SVs from lowest to highest coordinates on the reference
    2. Compute distances between adjacent SVs to get SV distance distribution
    3. Goodness-of-fit test with exponential distribution

    Returns
    -------
    DataFrame with two columns:
        chromosome name and the p-value from the test

    References
    ----------
    [1] Criteria for Inference of Chromothripsis in Cancer Genomes Korbel, Jan O. et al. Cell, Volume 152, Issue 6, 1226 - 1236
    """
    chr_groups = sv_df.group_by("chr")
    result: dict = {"chr": [], "clustering_significant": []}
    for chr, df in chr_groups:
        result["chr"].append(chr)
        if df.height > 2:
            df = df.sort("pos")
            distances: np.ndarray = (df["pos"][1:] - df["pos"][:-1]).to_numpy()
            mean = distances.mean()
            gof = goodness_of_fit(
                expon, known_params={"loc": 0, "scale": mean}, data=distances
            )
            if gof.fit_result.success and gof.pvalue < threshold:
                result["clustering_significant"].append(True)
                continue
        result["clustering_significant"].append(False)
    return pl.DataFrame(result)
