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

get_GSE207889 <- function(f, dir) {
  # [2026-07-02 Thu] ehh don't bother with this one. Very few hepatocytes
  library(Seurat)
  library(HGNChelper)
  library(openxlsx)
  obj <- CreateSeuratObject(counts = Read10X(here(dir, "GSM6822571")))
  obj <- NormalizeData(
    obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj, features = rownames(obj))
  obj <- RunPCA(obj, features = VariableFeatures(object = obj))
  obj <- FindNeighbors(obj, dims = 1:10)
  obj <- FindClusters(
    obj,
    resolution = 0.8,
    algorithm = 4,
    leiden_method = "igraph"
  )

  # Recommended preprocessing with seurat
  source(
    "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"
  )
  source(
    "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R"
  )
  gs_list <- gene_sets_prepare(
    "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx",
    "Liver"
  )
  score <- sctype_score(
    scRNAseqData = as.matrix(obj[["RNA"]]$scale.data),
    scaled = TRUE,
    gs = gs_list$gs_positive,
    gs2 = gs_list$gs_negative
  )
  annotation <- tibble(
    celltype = rownames(score)[max.col(t(score))],
    barcode = colnames(score)
  )
  # TODO: could wrap this up in another function. Figure out
  # what Seurat ScaleData does and see if scanpy can do something similar
}

## ** GSE154883
get_GSE154883 <- function(f, dir) {
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
  gse154883$write_h5ad(f)
}

## ** GSE264261

get_GSE264261 <- function(f, dir) {
  sc <- import("scanpy")
  adata <- sc$read_10x_h5(here(dir, "GSM8215174_filtered_feature_bc_matrix.h5"))
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
  adata$write_h5ad(f)
}

## ** GSE182604

get_GSE182604 <- function(f, dir) {
  df <- read_csv(here(dir, "GSE182604_PHH_scRNA_RawCount.csv.gz")) |>
    column_to_rownames(var = "...1")
  samples <- str_extract(colnames(df), "PHH.*") |> unique()

  adata <- lapply(samples, \(s) {
    cur <- df[, grepl(s, colnames(df))] |> t()
    obs <- data.frame(row.names = rownames(cur))
    obs$sample <- s
    ad$AnnData(X = cur, obs = obs)
  }) |>
    ad$concat(axis = "obs")
  adata$write_h5ad(f)
}


## ** Runner
extra_geo <- list(
  GSE154883 = get_GSE154883,
  GSE264261 = get_GSE264261,
  GSE182604 = get_GSE182604
)
for (gd in names(extra_geo)) {
  geo_file <- here(geo, gd, "adata.h5ad")
  if (!file.exists(geo_file)) {
    extra_geo[[gd]](geo_file, here(geo, gd))
  }
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
