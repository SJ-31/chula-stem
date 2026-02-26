suppressMessages({
  library(reticulate)
  library(tidyverse)
  library(paletteer)
  library(checkthat)
  library(glue)
  if (exists("snakemake")) {
    use_condaenv(snakemake@config$conda)
  } else {
    use_condaenv("stem-base")
  }
  ad <- import("anndata")
  sc <- import("scanpy")
  sc_utils <- import("chula_stem.sc_rnaseq")
  source(glue("{snakemake@config$r_src}/utils.R"))
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
  dc <- import("decoupler")

  adata <- ad$read_h5ad(snakemake@input[[1]])
  features <- read_lines(snakemake@input[[2]])
  adata <- adata[, rownames(adata$var) %in% features]
  prefix <- snakemake@params$prefix

  kws <- RCONFIG$kws

  bulked <- dc$pp$pseudobulk(adata, "sample", groups_col = NULL)
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
    RCONFIG$consensus_palette,
    prefix = prefix
  )
  tb <- lapply(seq(2, length(result)), \(k) {
    clst <- result[[k]]$consensusClass
    tibble(
      sample = names(clst),
      !!as.symbol(paste0(prefix, "_", k)) := clst
    )
  }) |>
    purrr::reduce(\(x, y) inner_join(x, y, by = join_by(sample)))

  write_csv(tb, snakemake@output[[1]])
  ggsave(snakemake@output[[2]], plot = plot, width = 25)
}

make_consensus_plot <- function(obj, prefix, palette = NULL) {
  plot <- ggplot(
    obj,
    aes(x = item, y = itemConsensus, fill = factor(cluster))
  ) +
    geom_bar(position = "fill", stat = "identity") +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("Sample") +
    ylab("Consensus score") +
    guides(fill = guide_legend(title = "Cluster assignment")) +
    facet_wrap(~k) +
    ggtitle(glue("Features: {prefix}"))
  if (!is.null(palette)) {
    plot <- plot + scale_fill_paletteer_d(palette)
  }
  plot
}

## * Visualizing DE genes & enriched pathways

# TODO: rename after the rule
prepare_tflink <- function() {
  id_map <- snakemake@config$gene_reference
  tflink_mitab <- RCONFIG$tflink_mitab
  tflink <- simplify_tflink_mitab(tflink_mitab)
  ncbi2hgnc <- setNames(id_map$hgnc_symbol, id_map$entrezgene_id)
  tflink$regulator <- ncbi2hgnc[tflink$regulator]
  tflink$target <- ncbi2hgnc[tflink$target]
  tflink <- filter(tflink, !(is.na(target) | is.na(regulator)))

  as_tbl_graph(tflink, directed = TRUE) |>
    activate(edges) |>
    mutate(ends = "last")
}

plot_de_regulatory <- function(G) {
  ggraph(G, "kk") +
    geom_edge_link(
      arrow = grid::arrow(
        ends = E(G)$ends,
        type = "closed",
        angle = 25,
        length = unit(0.1, "inches")
      ),
      end_cap = circle(0.75, "cm"),
      start_cap = circle(0.75, "cm")
    ) +
    geom_node_point(
      aes(fill = is_de, color = lfc, stroke = abs(lfc)),
      shape = 21,
      size = 5,
    ) +
    geom_node_label(aes(label = name), repel = TRUE) +
    theme(panel.background = element_rect(fill = "white")) +
    scale_color_paletteer_c("ggthemes::Green-Gold") +
    guides(
      fill = guide_legend("Is DE"),
      color = guide_legend("LFC")
    )
}

visualize_regulation <- function() {
  library(tidygraph)
  library(ggraph)

  tflink_g <- prepare_tflink()
}

## * Entry

if (exists("snakemake")) {
  globalenv()[[snakemake@rule]]()
}
