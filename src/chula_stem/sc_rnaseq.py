from collections import defaultdict
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Literal

import anndata as ad
import click
import numpy as np
import pandas as pd
import scanpy as sc
import skbio.stats.composition as comp
import sklearn.metrics as sm
from numba import jit
from scipy import sparse, stats

from chula_stem import utils as ut


@click.command()
@click.option(
    "-o", "--output", required=False, help="File to write counts to", default=""
)
@click.option(
    "-p",
    "--prefix",
    required=False,
    help="Prefix of files, e.g. <prefix>barcodes.tsv.gz",
    default="",
)
@click.option(
    "-d", "--dir", required=False, help="Directory containing files", default="."
)
def read_10x_mtx(output: str = "", prefix: str = "", dir: str = ".") -> ad.AnnData:
    adata: ad.AnnData = sc.read_10x_mtx(path=dir, prefix=prefix)
    if output and "h5ad" in output:
        adata.write(output)
    elif output:
        raise ValueError("Only writing to h5ad is supported")
    return adata


def make_marker_df(data: dict[str, list]) -> pd.DataFrame:
    # Create a sorted list of unique genes
    genes = sorted(set(gene for gene_list in data.values() for gene in gene_list))
    # Create a binary dataframe
    binary_df = pd.DataFrame(0, index=genes, columns=data.keys())
    # Fill the dataframe with 1 where the gene is expressed
    for cell_type, genes in data.items():
        binary_df.loc[genes, cell_type] = 1
    return binary_df


def cell_assign_wrapper(
    adata: ad.AnnData,
    cell_markers: pd.DataFrame | dict,
    model_path: Path = None,
    layer=None,
    size_factor_key="size_factors",
    train_kws: dict | None = None,
    **kws,
) -> dict:
    """Wrapper function for CellAssign from scvi

    Parameters
    ----------
    cell_markers : pd.DataFrame | dict
        Either a binary gene marker df of genes by cell types, or a dictionary mapping
        cell type names to lists of marker genes
    kws : kwargs
        kws passed to CellAssign

    Returns
    -------
    Dictionary of results for CellAssign
    Includes


    Notes
    -----

    """
    from scvi.external import CellAssign

    if isinstance(cell_markers, dict):
        cell_markers = (
            pd.DataFrame({"cell": cell_markers.keys(), "gene": cell_markers.values()})
            .explode("gene")
            .assign(value=1)
            .pivot_table(index="gene", columns="cell", values="value", fill_value=0)
            .astype(int)
        )
    filtered = adata[:, adata.var.index.isin(cell_markers.index)].copy()
    if model_path and model_path.exists():
        model = CellAssign.load(model_path)
    else:
        CellAssign.setup_anndata(filtered, size_factor_key=size_factor_key, layer=layer)
        model = CellAssign(filtered, cell_markers, **kws)
        model.train(**(train_kws or {}))
    predictions = model.predict()
    predictions.loc[:, "PREDICTION"] = predictions.idxmax(axis=1)
    predictions.index = adata.obs.index
    results = {"pred": predictions, "model": model}
    return results


def scib_normalize(
    adata,
    min_mean=0.1,
    log=True,
    precluster=True,
    clusters=(),
    cluster_method="louvain",
    sparsify=True,
):
    # A modified version of scib
    import anndata2ri
    import rpy2.rinterface_lib.callbacks
    import rpy2.rinterface_lib.embedded
    import rpy2.robjects as ro

    # Check for 0 count cells
    if np.any(adata.X.sum(axis=1) == 0):
        print("WARNING: found 0 count cells in the AnnData object")
        print("Filtering cells...")
        sc.pp.filter_cells(adata, min_counts=10, inplace=True)

    # Check for 0 count genes
    if np.any(adata.X.sum(axis=0) == 0):
        sc.pp.filter_genes(adata, min_counts=1, inplace=True)
        print("WARNING: found 0 count genes in the AnnData object")

    if sparsify:
        # massive speedup when working with sparse matrix
        if not sparse.issparse(adata.X):  # quick fix: HVG doesn't work on dense matrix
            adata.X = sparse.csr_matrix(adata.X)
    try:
        ro.r("library(scran)")
    except rpy2.rinterface_lib.embedded.RRuntimeError as ex:
        raise ImportError("scran library not installed")

    anndata2ri.activate()

    # keep raw counts
    adata.layers["counts"] = adata.X.copy()

    is_sparse = False
    x = adata.X.T
    # convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
    if sparse.issparse(x):
        is_sparse = True

        if x.nnz > np.power(2, 31) - 1:
            x = x.tocoo()
        else:
            x = x.tocsc()

    ro.globalenv["data_mat"] = x
    if len(clusters) > 0:
        precluster = False
    if precluster:
        # Preliminary clustering for differentiated normalisation
        adata_pp = adata.copy()
        sc.pp.normalize_total(adata_pp, counts_per_cell_after=1e6, min_counts=1)
        sc.pp.log1p(adata_pp)
        sc.pp.pca(adata_pp, n_comps=15, svd_solver="arpack")
        sc.pp.neighbors(adata_pp)
        if cluster_method == "louvain":
            sc.tl.louvain(adata_pp, key_added="groups", resolution=0.5)
        elif cluster_method == "leiden":
            sc.tl.leiden(adata_pp, key_added="groups", resolution=0.5)
        else:
            raise NotImplementedError(
                "Choose `cluster_method` from 'louvain', 'leiden'"
            )

        ro.globalenv["input_groups"] = adata_pp.obs["groups"]
        size_factors = ro.r(
            "sizeFactors("
            "   computeSumFactors("
            "       SingleCellExperiment(list(counts=data_mat)),"
            "       clusters = input_groups,"
            f"       min.mean = {min_mean}"
            "   )"
            ")"
        )

        del adata_pp
    elif len(clusters) > 0:
        ro.globalenv["clusters"] = clusters
        size_factors = ro.r(
            f"""
            sizeFactors(computeSumFactors(SingleCellExperiment(
            list(counts=data_mat)), min.mean = {min_mean}, clusters=clusters))
            """
        )
    else:
        size_factors = ro.r(
            f"""
            sizeFactors(computeSumFactors(SingleCellExperiment(
            list(counts=data_mat)), min.mean = {min_mean}))
            """
        )

    # modify adata
    adata.obs["size_factors"] = size_factors
    adata.X /= adata.obs["size_factors"].values[:, None]
    if log:
        print("Note! Performing log1p-transformation after normalization.")
        sc.pp.log1p(adata)
    else:
        print("No log-transformation performed after normalization.")

    if is_sparse:
        # convert to sparse, bc operation always converts to dense
        adata.X = sparse.csr_matrix(adata.X)

    adata.raw = adata  # Store the full data set in 'raw' as log-normalised data for statistical testing

    # Free memory in R
    ro.r("rm(list=ls())")
    ro.r("gc()")

    anndata2ri.deactivate()


def pca_to_leiden(
    adata,
    pca_kws=None,
    neighbor_kws=None,
    umap_kws=None,
    leiden_kws=None,
    do_umap: bool = False,
):
    """Wrapper function for doing the standard PCA to leiden clustering in scanpy"""
    if "X_pca" not in adata.obsm:
        ut.do_call(lambda **x: sc.pp.pca(adata, **x), pca_kws)
    ut.do_call(lambda **x: sc.pp.neighbors(adata, **x), neighbor_kws)
    if do_umap:
        ut.do_call(lambda **x: sc.tl.umap(adata, **x), umap_kws)
    ut.do_call(lambda **x: sc.tl.leiden(adata, **x), leiden_kws)


def annotate_marker(
    adata: ad.AnnData, marker_genes: list, gene_col: str = None, threshold: float = 0.0
) -> pd.DataFrame:
    """Mark cells where the marker gene expression is above the threshold (0 by default)"""
    expr: pd.DataFrame = (
        sc.get.obs_df(adata, keys=marker_genes, gene_symbols=gene_col) > threshold
    )
    for c in expr.columns:
        adata.obs[f"has_{c}"] = expr[c]


def annotate_adata_vars(
    adata: ad.AnnData,
    gene_id_col: str = "gene_id",
    savepath: Path | None = None,
    **kws,
) -> None:
    if not adata.var[gene_id_col].is_unique:
        raise ValueError(f"Can't pass a gene_id column '{gene_id_col}' with duplicates")
    if savepath and savepath.exists():
        metadata: pd.DataFrame = pd.read_csv(savepath)
    else:
        metadata = ut.get_ensembl_gene_data()
        metadata.to_csv(savepath)
    merged = (
        adata.var.reset_index(names="index")
        .merge(metadata, left_on=gene_id_col, right_on="ensembl_gene_id", how="left")
        .drop_duplicates("index")
        .set_index("index", drop=True)
    )
    adata.var = merged
    adata.var.loc[:, "mito"] = (adata.var["chromosome_name"] == "MT").fillna(False)


def distance_by_mads(
    adata: ad.AnnData,
    keys: Sequence | str | dict[str, str],
    group_keys: Sequence | str | None = None,
    inplace: bool = False,
) -> tuple[pd.DataFrame, dict] | None:
    """For each key, compute the distance of its values from the median in terms
        of its MAD i.e. (x - median(x)) / MAD(x)

    If `group_keys` is provided, by MAD and median are taken per group

    Notes
    ------
    Inspired by scuttle::isOutlier
    https://bioconductor.org/books/3.13/OSCA.basic/quality-control.html

    """
    if isinstance(keys, str):
        keys = [keys]
    result: dict = defaultdict(dict)

    def helper(cur: pd.DataFrame, group=None):
        tmp: dict = {}
        for key in keys:
            median = cur[key].median()
            mad = np.median(np.abs(cur[key] - median))
            if group:
                result[group][key] = mad.item()
            else:
                result[key] = mad.item()
            tmp[key] = (cur[key] - median) / mad
        return pd.DataFrame(tmp, index=cur.index)

    if group_keys is None:
        df = helper(adata.obs)
    else:
        dfs = []
        for group, idx in adata.obs.groupby(group_keys):
            cur = adata[idx, :].obs
            if cur.shape[0] > 0:
                cur_df = helper(cur, group=group)
                dfs.append(cur_df)
        df = pd.concat(dfs).loc[adata.obs.index, :]
    result = dict(result)
    if not inplace:
        return df, result
    adata.uns["mads"] = result
    adata.obsm["mads"] = df


# * Clustering utils


def neighbor_purity_from_distances(
    adata: ad.AnnData,
    cluster_key: str,
    distance_key: str = "distances",
    weight: bool = True,
    n_neighbors: int | None = None,
) -> np.ndarray:
    """Compute neighbor purity

    https://rdrr.io/bioc/bluster/man/neighborPurity.html
    """
    distances = adata.obsp[distance_key].astype(np.bool)
    clustering = adata.obs[cluster_key].astype(int)

    neighbor_assignments = distances.multiply(clustering.values).toarray()
    if n_neighbors is None and "neighbors" in adata.uns:
        n_neighbors = adata.uns["neighbors"]["params"]["n_neighbors"]
    else:
        raise ValueError("n_neighbors must be provided if not a key in adata.uns")

    neigbors_are_same = clustering.values.reshape(-1, 1) == neighbor_assignments
    n_in_cluster = neigbors_are_same.sum(axis=1)
    purity = n_in_cluster / n_neighbors
    if weight:
        weights = 1 / clustering.value_counts()
        weight_array = weights[clustering]
        return purity * weight_array.values
    return purity


def sweep_clustering(
    adata: ad.AnnData,
    cluster_fn: Callable,
    prefix: str,
    values: Sequence,
    distances: None | np.ndarray = None,
) -> None:
    """
    Perform a sweep over clustering parameter `param` across `values`,
    and assess the results with unsupervised clustering metrics

    Parameters
    ----------
    cluster_fn : Callable
        A function that takes two arguments: the anndata object and the parameter value
        and the name of the key in adata.obs to add the clustering results to
        e.g. lambda adata, res, key: sc.tl.leiden(adata, resolution = res, key_added = key)
        It must update the anndata object inplace

    prefix : str
        Prefix of clustering results key e.g. "leiden_res" for leiden_res1

    distances : np.ndarray | None
        Key in adata.obsp storing distances between cells e.g. `distances`
        after a call to sc.pp.neighbors

    TODO: would be nice if this were parallel
    """
    if distances is None:
        sc.pp.neighbors(adata)
        distances = adata.obsp["distances"]
    purity_tracker, sil_tracker = {}, {}
    for val in values:
        key_added = f"{prefix}{val}"
        cluster_fn(adata, val, key_added)
        result = adata.obs[key_added]

        purity = neighbor_purity_from_distances(adata, key_added, "distances", True)
        try:
            sil = sm.silhouette_samples(X=distances, labels=result)
            sil_mean = sil.mean()
        except ValueError:
            sil = sil_mean = np.nan
        metrics = {
            "silhouette_score": sil_mean,
            "mean_neighbor_purity": purity.mean(),
        }
        purity_tracker[key_added] = purity
        sil_tracker[key_added] = sil
        if key_added in adata.uns:
            adata.uns[key_added].update(metrics)
        else:
            adata.uns[key_added] = metrics
    sil_df = pd.DataFrame(sil_tracker).set_index(adata.obs_names)
    purity_df = pd.DataFrame(purity_tracker).set_index(adata.obs_names)
    adata.obsm[f"{prefix}_silhouette"] = sil_df
    adata.obsm[f"{prefix}_neighbor_purity"] = purity_df


def alr(
    data: ad.AnnData | pd.DataFrame,
    by: int | str,
    var_col: str | None = None,
    condition_col: str = "",
    zero_value: int = 1,
    gamma: float = 0,
) -> tuple[np.ndarray, int]:
    """Normalize counts in adata using ALR, with the counts of a specific
    gene `by` as the reference.

    Parameters
    ----------

    by : str | int
        name of gene to normalize by, or the index of the gene in adata.var
            if the name is provided, the index is looked up automatically.
    var_col : str | None
        column in adata.var containing the gene name. If not provided, adata.var_names
        is used
    """
    index: int | np.ndarray
    counts = data.X if isinstance(data, ad.AnnData) else data
    counts: np.ndarray = counts if isinstance(counts, np.ndarray) else counts.toarray()
    assert isinstance(counts, np.ndarray)
    counts[counts == 0] = zero_value
    if isinstance(by, str) and isinstance(data, pd.DataFrame):
        query = np.where(data.columns == by)
        index = query[0][0]
    elif isinstance(by, str) and var_col:
        query = np.where(data.var[var_col] == by)
        index = query[0][0]
        if len(query) > 1:
            raise ValueError("Key `by` is not unique!")
    elif isinstance(by, str):
        index = data.var.index.get_loc(by)
        if not isinstance(index, int):
            raise ValueError("Key `by` is not unique!")
    else:
        index = by
    return comp.alr(counts, index), index


@jit(nopython=True, cache=True)
def prop_helper(query_idx, all_idx, transformed):
    results = []
    v_cache: dict[int, float] = {}
    for i, q in enumerate(query_idx):
        q_vals = transformed[:, q]
        if q not in v_cache:
            v_cache[q] = q_vals.var()
        cur_var: float = v_cache[q]
        cur_result = []
        for other in all_idx:
            o_vals = transformed[:, other]
            if other not in v_cache:
                v_cache[other] = o_vals.var()
            o_var: float = v_cache[other]
            numer = (q_vals - o_vals).var()
            denom = cur_var + o_var
            if denom == 0:
                prop = np.nan
            else:
                prop: float = 1 - (numer / denom)
            cur_result.append(prop)
        results.append(np.array(cur_result))
    return results


def clr(
    x: ad.AnnData,
    features=None,
    feature_col: str = "GENEID",
    robust: bool = False,
    zero_value: int = 1,
) -> np.ndarray:
    "Centered log-ratio transform"
    counts: np.ndarray = x.X if isinstance(x.X, np.ndarray) else x.X.toarray()
    counts = counts + zero_value
    if features is None and not robust:
        return comp.clr(counts)
    if features is not None:
        subset = x.var[feature_col].isin(features)
    else:
        subset = np.ones(counts.shape[1]).astype(bool)
    if robust:
        gmean = [
            stats.gmean(counts[i, counts[i, :] != 0 & subset])
            for i in range(counts.shape[0])
        ]
        gmean = np.array(gmean)
    else:
        gmean = stats.gmean(counts[:, subset], axis=1, nan_policy="omit")
    clr = np.log(counts) - np.log(gmean.reshape(-1, 1))
    if robust:
        clr[clr == -np.inf] = 0
    return clr


def find_proportional_genes_rho(
    adata: ad.AnnData,
    query_genes: list[str],
    gene_col: str | None = None,
    lr_transform: Literal["clr", "alr"] = "clr",
    **kws,
) -> pd.DataFrame:
    if lr_transform == "clr":
        transformed: np.ndarray = clr(adata, **kws)
    elif lr_transform == "alr":
        transformed, removed_idx = alr(adata, **kws)
    tmp = adata.var_names if gene_col is None else pd.Index(adata.var[gene_col])

    if lr_transform == "alr":
        tmp = np.delete(tmp, removed_idx)

    all_genes: pd.Index = pd.Index(tmp)
    idx = list(range(len(all_genes)))
    results: list[np.ndarray] = prop_helper(
        [all_genes.get_loc(q) for q in query_genes],
        idx,
        transformed,
    )
    return pd.DataFrame(
        {q: results[i] for i, q in enumerate(query_genes)}, index=all_genes
    )
