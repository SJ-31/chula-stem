suppressMessages({
  library(reticulate)
  library(tidyverse)
  if (exists("snakemake")) {
    use_condaenv(snakemake@config$conda)
  } else {
    use_condaenv("stem-base")
  }
  ad <- import("anndata")
  sc <- import("scanpy")
  sc_utils <- import("chula_stem.sc_rnaseq")
})


if ("rlang" %in% utils::installed.packages()[, 1]) {
  options(error = rlang::entrace, rlang_backtrace_on_error = "full")
}

if (exists("snakemake")) {
  RCONFIG <- snakemake@config[[snakemake@rule]]
} else {
  RCONFIG <- list()
}

## * Cluster-fold similarity

write_cfs_results <- function(sce_list, sim_table, communities) {
  old_clusters <- lapply(
    sce_list,
    \(x) tibble(cell_id = colnames(x), cluster = colLabels(x))
  ) |>
    bind_rows()
  write_csv(old_clusters, snakemake@output$ind_clustering)

  sim_tb <- as_tibble(sim_table) |>
    mutate(
      comparison = apply(sim_table, 1, \(row) {
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

cfs_setup <- function(adata, subset_hvgs = FALSE, ...) {
  sc$pp$normalize_total(adata)
  sc$pp$log1p(adata)
  if (subset_hvgs) {
    sc$pp$highly_variable_genes(adata, subset = TRUE)
  }
  sc_utils$pca_to_leiden(adata, ...)
  sce <- AnnData2SCE(adata)
  assay(sce, "counts") <- assay(sce, "X")
  colLabels(sce) <- colData(sce)$leiden
  sce
}

do_cfs <- function(adata = NULL, kws = NULL, setup_kws = NULL) {
  library(ClusterFoldSimilarity)
  library(SingleCellExperiment)
  library(zellkonverter)

  if (is.null(adata)) {
    adata <- ad$read_h5ad(snakemake@input[[1]])
  }
  if (is.null(kws) && is.null(setup_kws)) {
    kws <- snakemake@config$do_cfs$cfs_kws
    setup_kws <- snakemake@config$do_cfs$setup_kws
  }
  adata <- adata[!adata$obs$sample %in% kws$cluster_exclude, ]
  samples <- unique(adata$obs$sample)
  sces <- lapply(
    samples,
    \(sample) {
      cur <- adata[adata$obs$sample == sample, ]$copy()
      do.call(\(...) cfs_setup(cur, ...), setup_kws)
    }
  ) |>
    `names<-`(samples)
  sim_table <- suppressMessages({
    do.call(
      \(...) clusterFoldSimilarity(sces, sampleNames = names(sces)),
      kws
    )
  })
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

## * Clustering

cluster_samples <- function() {
  library(ConsensusClusterPlus)

  adata <- ad$read_h5ad(snakemake@input[[1]])
  features <- read_lines(snakemake@input[[2]])
  adata <- adata[, rownames(adata$var) %in% features]
  prefix <- snakemake@params$prefix

  kws <- RCONFIG$kws

  bulked <- dc$pp$pseudobulk(tmp, "sample", groups_col = NULL)
  gene_mask <- sc$pp$filter_genes(bulked, min_counts = 50, inplace = FALSE)
  cell_mask <- sc$pp$filter_cells(bulked, min_counts = 50, inplace = FALSE)
  bulked <- bulked[cell_mask[[1]], gene_mask[[1]]]

  dge <- edgeR::DGEList(t(bulked$X))
  rownames(dge) <- rownames(bulked$var)
  colnames(dge) <- rownames(bulked$obs)
  dge <- edgeR::normLibSizes(dge, method = "TMM")
  expr <- edgeR::cpm(dge, log = TRUE)

  result <- do.call(\(...) ConsensusClusterPlus(d = expr, ...), kws)
  consensus <- calcICL(result)

  plot <- make_consensus_plot(
    consensus$itemConsensus,
    RCONFIG$consensus_palette
  )
  tb <- lapply(seq(2, length(result)), \(k) {
    clst <- result[[k]]$consensusClass
    tibble(sample = names(clst), !!as.symbol(k) := paste0(prefix, "_", clst))
  }) |>
    purrr::reduce(\(x, y) inner_join(x, y, by = join_by(sample)))

  write_csv(tb, snakemake@output[[1]])
  ggsave(plot, snakemake@output[[2]])
}

make_consensus_plot <- function(obj, palette = NULL) {
  plot <- ggplot(
    consensus$itemConsensus,
    aes(x = item, y = itemConsensus, fill = factor(cluster))
  ) +
    geom_bar(position = "fill", stat = "identity") +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("Sample") +
    ylab("Consensus score") +
    guides(fill = guide_legend(title = "Cluster assignment")) +
    facet_wrap(~k)
  if (!is.null(palette)) {
    plot <- plot + scale_fill_paletteer_d(palette)
  }
  plot
}

## * Entry

if (exists("snakemake")) {
  globalenv()[[snakemake@rule]]()
}
