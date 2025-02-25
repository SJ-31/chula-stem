#!/usr/bin/env ipython
import anndata as ad
import polars as pl
import polars.selectors as cs
import scanpy as sc
from chula_stem.sc_rnaseq import annotate_marker, pca_to_leiden, scib_normalize
from chula_stem.utils import read_existing
from pyhere import here

outdir = here("analyses", "output", "corneal", "tet_analysis")

MARKERS = ["TET1", "TET2", "TET3", "CD44", "MKI67"]

# Dataset: corneal
# Goal: Identify differentially expressed genes between MKI67 positive and negative
# clusters
# cells with high TET1 vs low TET1 expression
# Repeat for TET2,3

normalized_file = here(outdir, "corneal_normalized_no_cd44.h5ad")
normalized_wcd44_file = here(outdir, "corneal_normalized.h5ad")


def load(f):
    corneal_file = here(outdir.parent, "corneal.h5ad")
    corneal: ad.AnnData = ad.read_h5ad(corneal_file)
    corneal.obs.discard = ~corneal.obs.discard_reason.isin(
        ["kept", "high_subsets_mito_percent"]
    )
    corneal = corneal[~corneal.obs.discard, :]
    corneal = corneal[corneal.obs["scDblFinder.class"] != "doublet", :]
    annotate_marker(corneal, MARKERS)
    corneal.obs["MKI67_status"] = corneal.obs["has_MKI67"].replace(
        {True: "MKI67_positive", False: "MKI67_negative"}
    )
    pca_to_leiden(
        corneal,
        leiden_pars={"resolution": 1, "flavor": "igraph", "n_iterations": 2},
    )
    check_cd44 = (
        corneal.obs.groupby("leiden")
        .agg({"has_CD44": "sum", "leiden": "size"})
        .rename({"leiden": "cluster_size"}, axis="columns")
    )
    check_cd44["percent"] = check_cd44["has_CD44"] / check_cd44["cluster_size"]
    check_cd44.to_csv(here(outdir, "cd44_counts.csv"))

    scib_normalize(
        corneal,
        clusters=corneal.obs.leiden,
        precluster=False,
        cluster_method="",
    )

    no_cd44 = corneal[~corneal.obs["has_CD44"], :]

    for adata, outfile in zip(
        [corneal, no_cd44], [normalized_wcd44_file, normalized_file]
    ):
        sc.tl.rank_genes_groups(
            adata,
            groupby="MKI67_status",
            pts=True,
            rest=True,
            method="wilcoxon",
        )
        adata.write_h5ad(outfile)
    return ""


# Only compare fibroblasts
read_existing(normalized_file, load)

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
            locations = (self.names == gene).select(pl.all().arg_max())
            gene_values = []
            for indices, stats in zip(locations.iter_columns(), stat_df.iter_columns()):
                index = indices.first()
                if index:
                    gene_values.append(stats[index])
                else:
                    gene_values.append(None)
            results[gene] = gene_values
        return pl.DataFrame(results)


wcd44 = ad.read_h5ad(normalized_wcd44_file)
no_cd44 = ad.read_h5ad(normalized_file)
for adata in [no_cd44, wcd44]:
    ri = RankInterpreter(adata.uns["rank_genes_groups"])
    tet_percent = ri.gene_stats(MARKERS, "pts")
    tet_sig = ri.gene_stats(MARKERS, "pvals")
    tet_lfc = ri.gene_stats(MARKERS, "logfoldchanges")

    tet_sig.select(["group", MARKERS[0]]).filter(
        pl.any_horizontal(cs.by_dtype(pl.Float64) > 0.05)
    )

    # test = (ri.names == "MALAT1").select(pl.all().arg_max())


# Exclude non statistically significant groups

# High: must have at least 3rd quartile of lfc
# Low: must have at most 1st quartile of lfc

# #  --- CODE BLOCK ---
