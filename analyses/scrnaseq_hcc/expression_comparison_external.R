suppressMessages({
  library(here)
  library(reticulate)
  library(ggplot2)
  library(tidyverse)
  options(box.path = here("src"))
  box::use(R / plotting[dotplot])
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

## * Plotting

wanted_genes <- c(
  "BSG", # CD147
  "LIN28A",
  "LIN28B",
  "GPC3", # glypican 3
  "AFP" # Alpha-fetoprotein
)


combined <- ad$read_h5ad(combined_file)

norm_plot <- dotplot(
  obj = combined,
  var_names = wanted_genes,
  group_by = "sample",
  group_labels = list(
    "type" = "Redmonder::qMSOBu",
    "source" = "pals::glasbey"
  ),
  layer = "x_norm"
)
ggsave(
  here(data_dir, "dotplot_normalized.pdf"),
  plot = norm_plot
)

raw_plot <- dotplot(
  obj = combined,
  var_names = wanted_genes,
  group_by = "sample",
  group_labels = list(
    "type" = "Redmonder::qMSOBu",
    "source" = "pals::glasbey"
  )
)
ggsave(
  here(data_dir, "dotplot_raw.pdf"),
  plot = norm_plot
)
