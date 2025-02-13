#!/usr/bin/env ipython
import anndata as ad
import pandas as pd
import polars as pl
import polars.selectors as cs
import scanpy as sc
from chula_stem.utils import read_existing
from pyhere import here

outdir = here(
    "analyses",
    "output",
    "corneal",
)
from chula_stem.sc_rnaseq import annotate_marker, pca_to_leiden, scib_normalize

MARKERS = ["TET1", "TET2", "TET3"]

# Dataset: corneal
# Goal: Identify differentially expressed genes between
# cells with high TET1 vs low TET1 expression
# Repeat for TET2,3
corneal_normalized_file = here(outdir, "corneal_normalized.h5ad")


def load(f):
    corneal_file = here(outdir, "corneal.h5ad")
    corneal: ad.AnnData = ad.read_h5ad(corneal_file)
    corneal.obs.discard = ~corneal.obs.discard_reason.isin(
        ["kept", "high_subsets_mito_percent"]
    )
    corneal = corneal[~corneal.obs.discard, :]
    corneal = corneal[corneal.obs["scDblFinder.class"] != "doublet", :]
    pca_to_leiden(
        corneal, leiden_pars={"resolution": 1, "flavor": "igraph", "n_iterations": 2}
    )

    scib_normalize(
        corneal,
        clusters=corneal.obs.leiden,
        precluster=False,
        cluster_method="",
    )
    # <2025-02-13 Thu> Operating on leiden clusters for now, but you
    # can switch to cell types if you have good categorization of them
    sc.tl.rank_genes_groups(
        corneal, groupby="leiden", pts=True, rest=True, method="wilcoxon"
    )
    annotate_marker(corneal, MARKERS)
    corneal.write_h5ad(f)


adata = read_existing(corneal_normalized_file, load, ad.read_h5ad)

# Plot marker expression
# fig = sc.pl.umap(adata, color = )

# #  --- CODE BLOCK ---


class RankInterpreter:
    """A class with methods to ease data access and interpretation of the results
    produced by scanpy.tl.rank_genes_groups
    """

    def __init__(self, results: dict) -> None:
        for k, v in results.items():
            if k != "params":
                self.__setattr__(k, pl.DataFrame(v))
        self.groups = self.names.columns
        # attrs include 'params', 'pts' (optional), 'pts_rest' (optional),
        # 'names',
        # 'scores',
        # 'pvals',
        # 'pvals_adj',
        # 'logfoldchanges'

    def gene_stats(self, gene_list, stat: str) -> pl.DataFrame:
        """Get a dataframe of group x genes
        showing the values of each gene in `gene_list` for the given statistic.
        """
        stat_df: pl.DataFrame = self.__getattribute__(stat)
        results: dict = {"group": self.groups}
        for gene in gene_list:
            mask: pl.DataFrame = self.names == gene
            index_expr = [
                pl.when(pl.col(c)).then(pl.int_range(pl.len())).otherwise(None).alias(c)
                for c in mask.columns
            ]
            with_index = mask.with_columns(index_expr).min()
            gene_values = []
            for indices, stats in zip(
                with_index.iter_columns(), stat_df.iter_columns()
            ):
                index = indices.first()
                if index:
                    gene_values.append(stats[index])
                else:
                    gene_values.append(None)
            results[gene] = gene_values
        return pl.DataFrame(results)


ri = RankInterpreter(adata.uns["rank_genes_groups"])
tet_percent = ri.gene_stats(MARKERS, "pts")
tet_sig = ri.gene_stats(MARKERS, "pvals")
tet_lfc = ri.gene_stats(MARKERS, "logfoldchanges")

tet_sig.select(["group", MARKERS[0]]).filter(
    pl.any_horizontal(cs.by_dtype(pl.Float64) > 0.05)
)
# Exclude non statistically significant groups

# High: must have at least 3rd quartile of lfc
# Low: must have at most 1st quartile of lfc

# #  --- CODE BLOCK ---
