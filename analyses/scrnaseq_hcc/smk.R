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
  ENV <- snakemake@config
} else {
  RCONFIG <- list()
  ENV <- list()
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
    read_csv(ENV$gene_reference),
    ENV$tflink_mitab
  )
  fi_g <- prepare_reactome_fi(ENV$reactome_fi)
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
        ggsave_graph_dynamic(
          tflink_to_plot,
          plot_de_graph(tflink_to_plot, palette_c = RCONFIG$palette_c),
          glue("{outdir}/{prefix}_tflink_{key}.pdf")
        )
      }

      fi_to_plot <- annotate_graph_de(
        fi_g,
        cur_cluster,
        min_interesting = 2,
        kept_nodes = RCONFIG$always_include
      )
      if (!is.null(fi_to_plot)) {
        ggsave_graph_dynamic(
          fi_to_plot,
          plot_de_graph(fi_to_plot, palette_c = RCONFIG$palette_c),
          glue("{outdir}/{prefix}_fi_{key}.pdf")
        )
      }
    }
  }
}


## ** Enrichment pathways

plot_go_graph <- function(tb, go, output_file) {
  comps <- to_go_graph_components(
    tb,
    go,
    min_enriched = RCONFIG$min_enriched,
    simplify_threshold = RCONFIG$simplify_threshold,
    context_threshold = RCONFIG$context_threshold,
    min_ns_dist = RCONFIG$min_ns_dist
  ) |>
    lapply(\(g) {
      activate(g, nodes) |>
        mutate(
          name = ifelse(is.na(name), name, str_wrap(glue("{name}\n{term}"), 30))
        )
    })
  tmp <- tempdir()
  files <- lapply(
    seq_along(comps),
    \(i) {
      plot <- plot_enriched_graph(
        comps[[i]],
        palette_c = RCONFIG$palette_c %||% "ggthemes::Classic Red-Green Light",
        palette_d = RCONFIG$palette_d %||% "LaCroixColoR::CeriseLimon"
      )
      name <- glue("{tmp}/{i}.pdf")
      ggsave_graph_dynamic(comps[[i]], plot, name)
      list(name) |> `names<-`(length(comps[[i]]))
    },
  ) |>
    list_c()
  files <- files[order(names(files))]
  join_pdfs(files, output_file)
}

plot_reactome_graph <- function(tb, rg, output_file) {
  comps <- rg |>
    activate(nodes) |>
    left_join(tb, by = join_y(name)) |>
    mutate(enriched = ifelse(enriched, TRUE, FALSE)) |>
    keep_interesting_comps(rg, "enriched")
  tmp <- tempdir()
  files <- lapply(
    seq_along(comps),
    \(i) {
      cur <- keep_interesting_leaves(comps[[i]], "enriched") |>
        activate(nodes) |>
        mutate(
          name = str_wrap(name, width = max_length)
        )
      plot <- plot_enriched_graph(
        cur,
        palette_c = RCONFIG$palette_d %||% "ggthemes::Classic Red-Green Light",
        palette_d = RCONFIG$palette_d %||% "LaCroixColoR::CeriseLimon"
      )
      name <- glue("{tmp}/{i}.pdf")
      ggsave_graph_dynamic(cur, plot, name)
      list(name) |> `names<-`(length(cur))
    }
  ) |>
    list_c()
  files <- files[order(names(files))]
  join_pdfs(files, output_file)
}

visualize_enrichment_inner <- function(outdir) {
  GO <- read_graph(ENV$go_graph, "gml") |>
    as_tbl_graph() |>
    select(-id) |> # BUG: fix the name weirdness
    rename(distance_to_ns = "distancetons")
  go_term2id <- read_csv(ENV$go_term2id) |> distinct(term, .keep_all = TRUE)
  RG <- prepare_reactome_pwy(
    pathway_file = ENV$reactome_pathways,
    relation_file = ENV$reactome_relations
  )
  reactome_ids2name <- read_tsv(
    ENV$reactome_pathways,
    col_names = c("id", "name", "species")
  ) |>
    select(id, name) |>
    deframe()

  # Combine them, loop over, and get results
  inputs <- list(
    gprofiler_clusters = read_csv(snakemake@input$gprofiler_clusters),
    gprofiler_samples = read_csv(snakemake@input$gprofiler_samples),
    decoupler_clusters = read_csv(snakemake@input$clusters_gs)
  )
  gprofiler_renaming <- c(
    name = "native",
    pvalue = "p_value",
    stat = "precision",
    enriched = "significant"
  )

  for (n in names(inputs)) {
    tb <- inputs[[n]]
    if (str_starts(n, "gprofiler")) {
      go_tb <- tb |>
        filter(str_starts(native, "GO:")) |>
        dplyr::select(native, precision, p_value) |>
        rename(gprofiler_renaming)

      reactome_tb <- tb |>
        filter(str_starts(native, "REAC:")) |>
        mutate(
          native = str_remove(native, "REAC:"),
          native = reactome_ids2name[native]
        ) |>
        dplyr::filter(!is.na(native)) |>
        dplyr::select(native, precision, p_value) |>
        rename(gprofiler_renaming)
    } else {
      go_tb <- tb |>
        filter(str_starts(name, "GO_..:")) |>
        mutate(name = str_remove(name, "GO_..:"), enriched = padj <= 0.05) |>
        inner_join(go_term2id, by = join_by(x$name == y$term))
      reactome_tb <- tb |>
        filter(str_starts(name, "Reactome:")) |>
        mutate(name = str_remove(name, "^Reactome:"), enriched = padj <= 0.05)
    }
    group_key <- ifelse(str_ends(n, "clusters"), "clustering", "analysis_group")

    groupings <- tb |>
      group_by(!!as.symbol(group_key), contrast) |>
      summarise() |>
      deframe()

    for (i in seq_along(groupings)) {
      gk <- names(groupings[i])
      cur_contrast <- groupings[i]

      go_cur <- go_tb |>
        filter(contrast == cur_contrast, !!as.symbol(group_key) == gk)
      reactome_cur <- reactome_tb |>
        filter(contrast == cur_contrast, !!as.symbol(group_key) == gk)

      prefix <- glue("{gk} {cur_contrast}") |>
        str_replace_all(" ", "_")
      suffix <- ifelse(str_ends(n, "clusters"), "clusters", "samples")

      reactome_out <- glue("{outdir}/{prefix}_reactome_{suffix}.pdf")
      go_out <- glue("{outdir}/{prefix}_go_{suffix}.pdf")

      plot_go_graph(tb = go_tb, go = GO, output_file = go_out)
      plot_reactome_graph(tb = reactome_tb, rg = RG, output_file = reactome_out)
    }
  }
}

visualize_enrichment <- function() {
  outdir <- snakemake@params$outdir
  dir.create(outdir, showWarnings = FALSE)
  tryCatch(
    expr = visualize_enrichment_inner(outdir),
    error = \(.) {
      unlink(outdir, recursive = TRUE)
    }
  )
}

## * Entry

if (exists("snakemake")) {
  globalenv()[[snakemake@rule]]()
}
