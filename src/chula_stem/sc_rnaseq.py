import anndata as ad
import click
import pandas as pd
import scanpy as sc
from scvi.external import CellAssign


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
