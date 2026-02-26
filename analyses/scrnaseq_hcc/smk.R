suppressMessages({
  library(reticulate)
  library(tidyverse)
  library(paletteer)
  library(checkmate)
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
# TODO: you should save this somewhere else

annotate_graph_de <- function(G, tb, min_interesting = 1, kept_nodes = NULL) {
  kept <- activate(G, nodes) |>
    left_join(tb, by = join_by(x$name == y$gene)) |>
    mutate(is_de = !is.na(lfc)) |>
    keep_interesting_comps("is_de", min_interesting, kept_nodes = kept_nodes)
  if (length(kept) > 0) {
    purrr::reduce(kept, \(x, y) bind_graphs(x, y))
  } else {
    NULL
  }
}

visualize_de_genes <- function() {
  source(glue("{snakemake@config$r_src}/plotting.R"))
  dir.create(snakemake@params$outdir)
  library(tidygraph)
  library(igraph)
  library(ggraph)

  tflink_g <- prepare_tflink(
    read_csv(snakemake@config$gene_reference),
    RCONFIG$tflink_mitab
  )
  fi_g <- prepare_reactome_fi(RCONFIG$reactome_fi)
  clusters_de <- local({
    scvi <- read_csv(snakemake@input$cluster_level_scvi) |>
      filter(proba_de > RCONFIG$scvi_min_prob) |>
      rename(lfc = "lfc_median")
    edger <- read_csv(snakemake@input$cluster_level_edger) |>
      rename(lfc = "logFC")
    bind_rows(edger, scvi) |>
      select(gene, lfc, contrast, clustering) |>
      distinct(gene, contrast, clustering, .keep_all = TRUE)
  })
  samples_de <- read_csv(snakemake@input$sample_level)

  outdir <- snakemake@params$outdir
  key2tb <- list(clusters = clusters_de, samples = samples_de)
  for (key in names(key2tb)) {
    group_key <- ifelse(key == "clusters", "clustering", "analysis_group")
    tb <- key2tb[[key]]

    groupings <- key2tb[[key]] |>
      group_by(!!as.symbol(group_key), contrast) |>
      summarise() |>
      deframe()

    for (i in seq_along(groupings)) {
      cur_contrast <- groupings[i]
      gk <- names(groupings[i])

      prefix <- glue("{gk} {cur_contrast}") |>
        str_replace_all(" ", "_")
      cur_cluster <- tb |>
        filter(
          !!as.symbol(group_key) == gk & contrast == cur_contrast
        )

      tflink_to_plot <- annotate_graph_de(
        tflink_g,
        cur_cluster,
        min_interesting = 1,
        kept_nodes = RCONFIG$always_include
      )
      if (!is.null(tflink_to_plot)) {
        ggsave(
          glue("{outdir}/{prefix}_tflink_{key}.pdf"),
          plot = plot_de_graph(tflink_to_plot),
          height = 15,
          width = 15
        )
      }

      fi_to_plot <- annotate_graph_de(
        fi_g,
        cur_cluster,
        min_interesting = 2,
        kept_nodes = RCONFIG$always_include
      )
      if (!is.null(fi_to_plot)) {
        ggsave(
          glue("{outdir}/{prefix}_fi_{key}.pdf"),
          plot = plot_de_graph(fi_to_plot),
          height = 15,
          width = 15
        )
      }
    }
  }
}

## * Entry

if (exists("snakemake")) {
  globalenv()[[snakemake@rule]]()
}
