suppressMessages({
  library(here)
  library(reticulate)
  library(ggplot2)
  library(tidyverse)
  options(box.path = here("src"))
  suppressWarnings(box::use(R / plotting[dotplot]))
  use_condaenv("stem-base")
  ad <- import("anndata")
})

data <- here("analyses", "data_all")
data_dir <- here(
  data,
  "output",
  "HCC",
  "SCRNASEQ",
  "cohort",
  "compare_external"
)
geo <- here(data, "public_data", "GEO")
gene_reference <- here("analyses/data/ensembl_gene_data.csv")
combined_file <- here(data_dir, "combined.h5ad")

## * Collect misc. data

## ** GSE154883
gse154883_file <- here(geo, "GSE154883", "adata.h5ad")
if (!file.exists(gse154883_file)) {
  library(Seurat)
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

## ** GSE264261

gse264261_file <- here(geo, "GSE264261", "adata.h5ad")
if (!exists()) {
  sc <- import("scanpy")
  adata <- sc$read_10x_h5(here(
    geo,
    "GSE264261",
    "GSM8215174_filtered_feature_bc_matrix.h5"
  ))
  obs <- read_tsv(here(geo, "GSE264261", "GSM8215174_clusters.tsv.gz")) |>
    column_to_rownames(var = "barcode")
  adata$obs <- merge(adata$obs, obs, by = "row.names") |>
    column_to_rownames("Row.names")
  adata <- adata[adata$obs$status != "doublet", ]
  adata$obs$sample <- ifelse(
    adata$obs$assignment == 0,
    "GSE264261_cluster0",
    "GSE264261_cluster1"
  )
  adata$obs <- select(adata$obs, sample)
  adata$write_h5ad(gse264261_file)
}

## * Plotting

wanted_genes <- list(
  main = c(
    "BSG", # CD147
    "LIN28A",
    "LIN28B",
    "GPC3", # glypican 3
    "AFP"
  ), # Alpha-fetoprotein
  edgeR_DE = c(
    # Genes DE with edgeR
    "SDF2L1",
    "NUDT14",
    "SCAND1",
    "MRPL55",
    "RPL28",
    "C4orf48",
    "HCFC1R1",
    "ABHD17A",
    "GRINA",
    "MEA1"
  ),
  scVI_DE = c(
    "GAL",
    "TAGLN",
    "MYL9",
    "RFLNB"
  )
)

group_labels <- list(
  "type" = "Redmonder::qMSOBu",
  "source" = "pals::glasbey"
)

combined <- ad$read_h5ad(combined_file)

save_dot <- function(plot, fname) {
  ggsave(fname, width = 12, height = 9, plot = plot)
}

dotplot(
  obj = combined,
  var_names = wanted_genes,
  group_by = "sample",
  group_labels = group_labels,
  layer = "x_norm",
  sort = "type"
) |>
  save_dot(here(data_dir, "dotplot_normalized.pdf"))

dotplot(
  obj = combined,
  var_names = wanted_genes,
  group_by = "sample",
  group_labels = group_labels,
  sort = "type"
) |>
  save_dot(here(data_dir, "dotplot_raw.pdf"))
