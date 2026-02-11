library(reticulate)
library(tidyverse)
library(zellkonverter)
library(SingleCellExperiment)
use_condaenv("stem-base")
ad <- import("anndata")


## * Cluster-fold similarity

write_cfs_results <- function(sce_list, sim_table, communities) {
  old_clusters <- lapply(
    sce_list,
    \(x) tibble(cell_id = colnames(x), cluster = colLabels(x))
  ) |>
    bind_rows()
  write_csv(old_clusters, snakemake@output$independent_clustering)

  sim_tb <- as_tibble(sim_table) |>
    mutate(
      comparison = apply(res, 1, \(row) {
        paste0(
          sort(c(
            paste0(row["datasetL"], row["clusterL"]),
            paste0(row["datasetR"], row["clusterR"])
          )),
          collapse = ";"
        )
      })
    ) |>
    distinct(comparison, .keep_all = TRUE) |>
    select(-comparison) |>
    arrange(desc(similarityValue))

  write_csv(sim_tb, snakemake@output$similarity)
  write_csv(communities, snakemake@output$communities)
}

cfs_setup <- function(adata, subset_hvgs = FALSE) {
  sc$pp$normalize_total(adata)
  sc$pp$log1p(adata)
  if (subset_hvgs) {
    sc$pp$highly_variable_genes(adata, subset = TRUE)
  }
  sc_utils$pca_to_leiden(adata)
  sce <- AnnData2SCE(adata)
  assay(sce, "counts") <- assay(sce, "X")
  colLabels(sce) <- colData(sce)$leiden
  sce
}

do_cfs <- function(adata, kws, setup_kws) {
  library(ClusterFoldSimilarity)

  if (is.null(adata)) {
    adata <- ad$read_h5ad()
  }
  if (is.null(kws)) {
    kws <- snakemake@config$do_cfs$cfs_kws
    setup_kws <- snakemake@config$do_cfs$setup_kws
  }
  samples <- unique(adata$obs$sample)
  sces <- lapply(
    samples,
    \(sample) {
      cur <- adata[adata$obs$sample == sample, ]$copy()
      do.call(\(...) cfs_setup(cur, ...), setup_kws)
    }
  ) |>
    `names<-`(samples)
  sim_table <- do.call(
    \(...) clusterFoldSimilarity(sces, sampleNames = names(sces)),
    kws
  )
  communities <- findCommunitiesSimmilarity(sim_table)

  if (exists("snakemake")) {
    write_cfs_results(
      sces,
      sim_table = sim_table,
      communities = communities
    )
  }
  list(similarities = sim_table, communities = communities)
}


## * Entry

if (exists("snakemake")) {
  globalenv()[[snakemake@rule]]()
}
