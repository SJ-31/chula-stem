library(here)
library(glue)
library(tidyverse)
library(ConsensusClusterPlus)
library(edgeR)
library(reticulate)
use_condaenv("stem-base")
source(here("src", "R", "de_analysis.R"))

yte <- import("yte")
env <- yte$process_yaml(paste0(
  read_lines(here("analyses", "cluster_reproducibility", "env.yaml")),
  collapse = "\n"
))
seed <- env$seed

ad <- import("anndata")
sc <- import("scanpy")

workdir <- do.call(here, as.list(env$workdir))
sim <- do.call(here, as.list(env$outs$simulated))
outdir <- do.call(here, as.list(env$outs$clustered))

dir.create(outdir, recursive = TRUE)

reals <- map_chr(env$datasets$real, \(x) do.call(here, as.list(x)))

adatas <- c(list.files(sim, full.names = TRUE), reals)

cluster_call_de <- function(file) {
  library(ConsensusClusterPlus)

  adata <- ad$read_h5ad(file)
  sc$pp$filter_cells(adata, min_genes = 50)
  sc$pp$filter_cells(adata, min_counts = 100)
  sc$pp$filter_genes(adata, min_cells = 20)

  if (adata$shape[[1]] == 0) {
    warning(glue("File {basename(file)} has no count data after filtering"))
    return()
  } else {
    message(glue("Proceeding {basename(file)} to DE analysis"))
  }

  x <- adata$X
  if (!is.matrix(x)) {
    x <- as.matrix(x)
    x[is.na(x)] <- 0
  }
  dge <- edgeR::DGEList(t(x) + 1)

  dge$genes$gene <- rownames(adata$var)

  dge <- edgeR::normLibSizes(dge, method = "TMM")
  expr <- edgeR::cpm(dge, log = TRUE)

  result <- ConsensusClusterPlus(d = expr, maxK = 7)
  ic_vals <- calcICL(result)
  best <- as_tibble(ic_vals$clusterConsensus) |>
    group_by(k) |>
    summarise(avg = mean(clusterConsensus, na.rm = TRUE)) |>
    slice_max(avg) |>
    pluck("k")

  tb <- lapply(seq(2, length(result)), \(k) {
    clst <- result[[k]]$consensusClass
    tibble(
      sample = names(clst),
      !!as.symbol(paste0("cc_", k)) := clst
    )
  }) |>
    purrr::reduce(\(x, y) inner_join(x, y, by = join_by(sample))) |>
    column_to_rownames("sample")

  best_cc <- tb[[glue("cc_{best}")]]

  dge$samples$best_cc <- as.character(best_cc)
  dge$samples$sample <- rownames(adata$obs)
  de_result <- edgeR_glm_wrapper(dge, group = "best_cc", id_col = "gene")

  adata$obs <- tb
  adata$obs$best_cc <- best_cc
  adata$uns$clusterConsensus <- as.data.frame(ic_vals$clusterConsensus)
  adata$uns$itemConsensus <- as.data.frame(ic_vals$itemConsensus)
  adata$uns$top_de <- de_result$top
  adata
}

for (file in adatas) {
  out <- here(outdir, basename(file))
  if (!file.exists(out)) {
    result <- cluster_call_de(file)
    if (!is.null(result)) {
      result$write_h5ad(out)
    }
  }
}
