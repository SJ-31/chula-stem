import anndata as ad
import click
import scanpy as sc


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
