library(here)
library(reticulate)
library(tidyverse)
library(Seurat)

use_condaenv("stem-base")

ad <- import("anndata")

data <- here("analyses", "data_all")
outdir <- here("analyses", "output", "scrnaseq_hcc")
geo <- here(data, "public_data", "GEO")
gene_reference <- here("analyses/data/ensembl_gene_data.csv")

gse154883_file <- here(geo, "GSE154883", "adata.h5ad")
if (!file.exists(gse154883_file)) {
  load(here(
    geo,
    "GSE154883",
    "GSE154883_ARPKD.hashtag.HTOdemux.Robj"
  ))
  imported <- UpdateSeuratObject(arpkd.hashtag)
  imported <- imported[, imported[[]]$HTO_classification.global == "Singlet"]
  imported <- imported[, grepl("Ctrl", imported[[]]$HTO_maxID)]

  var <- tibble(hgnc_symbol = Features(imported)) |>
    left_join(
      distinct(
        select(read_csv(gene_reference), ensembl_gene_id, hgnc_symbol),
        hgnc_symbol,
        .keep_all = TRUE
      ),
      by = join_by(hgnc_symbol)
    ) |>
    rename(gene_ids = "ensembl_gene_id") |>
    column_to_rownames("hgnc_symbol")

  gse154883 <- ad$AnnData(
    X = t(LayerData(imported, "counts")),
    var = var,
    obs = data.frame(
      source = "GSE154883",
      sample = imported[[]]$HTO_maxID,
      row.names = colnames(imported)
    )
  )
  gse154883$write_h5ad(gse154883_file)
}
