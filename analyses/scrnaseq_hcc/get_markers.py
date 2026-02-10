#!/usr/bin/env ipython


import polars as pl
from pyhere import here

too_predict = here().parent / "too-predict"

id_map = pl.read_csv(
    too_predict / "/data/reference/Homo_sapiens.GRCh38.113.gene_id_mapping.tsv",
    separator="\t",
).rename({"gene_name": "hgnc_symbol"})
cell_markers_old = pl.read_csv(
    too_predict / "data/reference/cell_markers_custom_meta.tsv",
    separator="\t",
)
gene_sets_old = pl.read_csv(
    too_predict / "data/reference/gene_sets_custom_meta.tsv",
    separator="\t",
)

# Normal may also have cholangiocyte

cell_markers = {
    "t_helper_17_cell": [],
    "gamma_delta_t_cell": [],
    "m1_macrophage": [],
    "t_helper_1_cell": [],
    "m2_macrophage": [],
    "exhausted_t_cell": [],
    "stromal_cell": [],
    "natural_killer_cell": [],
    "cancer_cell": [],
    "neutrophil": [],
    "exhausted_cd8+_t_cell": [],
    "macrophage": [],
    "monocyte": [],
    "hepatocyte": [],
    # TODO: find for cholangiocytes
}
for names, group in (
    cell_markers_old.join(id_map, left_on="ensembl", right_on="gene_id")
    .filter(pl.col("hgnc_symbol").is_not_null())
    .group_by("cell_type")
):
    if names[0] in cell_markers:
        cell_markers[names[0]].extend(list(set(group["hgnc_symbol"])))


# TODO: rename the genes in gene sets
