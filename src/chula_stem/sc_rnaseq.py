import anndata
import anndata as ad
import click
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scvi.external import CellAssign

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
    marker_df: pd.DataFrame,
    model_path: str,
    type_key: str = "cell_type",
    count_layer="counts",
    size_factor_key="size_factors",
) -> pd.DataFrame:
    filtered = adata[:, adata.var.index.isin(marker_df.index)].copy()
    try:
        model = CellAssign.load(model_path)
    except:
        CellAssign.setup_anndata(
            filtered, size_factor_key=size_factor_key, layer=count_layer
        )
        model = CellAssign(filtered, marker_df)
    model.train()
    model.save(model_path, save_anndata=True, overwrite=True)
    predictions = model.predict()
    adata.obs[type_key] = predictions.idxmax(axis=1).values
    return predictions


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
        sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6, min_counts=1)
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
    adata, pca_pars=None, neighbor_pars=None, umap_pars=None, leiden_pars=None
):
    """Wrapper function for doing the standard PCA to leiden clustering in scanpy"""
    if "X_pca" not in adata.obsm:
        ut.do_call(lambda **x: sc.pp.pca(adata, **x), pca_pars)
    ut.do_call(lambda **x: sc.pp.neighbors(adata, **x), neighbor_pars)
    ut.do_call(lambda **x: sc.tl.umap(adata, **x), umap_pars)
    ut.do_call(lambda **x: sc.tl.leiden(adata, **x), leiden_pars)


def annotate_marker(
    adata: ad.AnnData, marker_genes: list, gene_col: str = None
) -> pd.DataFrame:
    """Mark cells where the marker gene expression is nonzero"""
    expr: pd.DataFrame = (
        sc.get.obs_df(adata, keys=marker_genes, gene_symbols=gene_col) > 0
    )
    for c in expr.columns:
        adata.obs[f"has_{c}"] = expr[c]
    return expr.agg("sum")
