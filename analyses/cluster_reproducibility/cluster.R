library(here)
library(glue)
library(tidyverse)
library(ConsensusClusterPlus)
library(edgeR)
library(SPsimSeq)
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

workdir <- do.call(here, as.list(env$workdir))
outdir <- do.call(here, as.list(env$outs$clustered))
datasets <- here(outdir, "datasets")

dir.create(outdir)

adatas <- c(list.files(datasets), env$datasets$real)

cluster_call_de <- function(file) {
  library(ConsensusClusterPlus)

  adata <- ad$read_h5ad(file)

  dge <- edgeR::DGEList(adata$X)
  rownames(dge) <- rownames(adata$var)
  dge$genes$gene <- rownames(adata$var)
  colnames(dge) <- rownames(adata$obs)

  dge <- edgeR::normLibSizes(dge, method = "TMM")
  expr <- edgeR::cpm(dge, log = TRUE)

  result <- ConsensusClusterPlus(d = expr)
  ic_vals <- calcICL(result)
  best <- as_tibble(ic_vals$clusterConsensus) |>
    group_by(k) |>
    summarise(avg = mean(ClusterConsensus, na.rm = TRUE)) |>
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

  dge$samples$best_cc <- best_cc
  de_result <- edgeR_glm_wrapper(dge, group = "best_cc", id_col = "gene")

  adata$obs <- merge(adata$obs, tb, by = "row.names")
  adata$obs$best_cc <- best - cc
  adata$uns$clusterConsensus <- as.data.frame(ic_vals$clusterConsensus)
  adata$uns$itemConsensus <- as.data.frame(ic_vals$itemConsensus)
  adata$uns$top_de <- de_result$top
  adata
}

for (file in adatas) {
  result <- cluster_call_de(file)
}
