R_SRC <- Sys.getenv("R_SRC")
source(paste0(R_SRC, "/", "utils.R"))
library(scDblFinder)
library(scater)
library(edgeR)
library(tidyverse)
library(glue)
library(zellkonverter)

#' Compute Gini impurity of clustered SingleCellExperiment object
#'  as described in scBubbleTree package
#'
#' @param clusters either a vector of cluster assignments for each cell or a column
#'    of colData(sce) containing those assignments
#' @param label_col a column of colData(sce) to compute Gini impurity on
sce_gini <- function(sce, clusters, label_col) {
  df <- colData(sce)
  if (length(clusters) == ncol(df)) {
    df$temp.clusters <- clusters
    clusters <- "temp.clusters"
  }
  if (!clusters %in% colnames(df)) {
    stop("`clusters` must be a valid column in sce's colData!")
  }
  ginis <- list()
  n <- nrow(df)
  whole <- lapply(unique(df[[clusters]]), \(x) {
    cur <- df[df[[clusters]] == x, ]
    n_i <- nrow(cur)
    freq <- table(cur[[label_col]]) / n_i
    gini <- sum(freq * (1 - freq))
    ginis[[x]] <<- gini
    gini * (n_i / n) # Apply weight
  }) |>
    unlist() |>
    sum()
  list(cluster_gini = ginis, whole_gini = whole)
}


add_feature_info <- function(sce, db, keytype = "GENEID") {
  if (is.character(db)) {
    db <- EnsDb(db)
  }

  anno_cols <- c("GENENAME", "GENEID", "GENEBIOTYPE", "SEQNAME")
  anno <- AnnotationDbi::select(db,
    keys = rownames(sce), keytype = keytype,
    columns = anno_cols
  ) |> as_tibble()
  tmp <- tibble(join = rownames(sce))

  anno$join <- anno[[keytype]]

  rowData(sce) <- left_join(tmp, anno, by = join_by(join)) |>
    dplyr::distinct(join, .keep_all = TRUE) |>
    dplyr::select(-join)

  n_missing <- is.na(rowData(sce)$GENEID) |> sum()
  message(glue("Missing genes: {n_missing}"))
  rowData(sce)$is_mito <- (rowData(sce)$SEQNAME == "MT") %>% replace_na(FALSE)
  sce
}

# Requires a SingleCellExperiment object
# Use zellkonverter as entry into SCverse

#' Read count data into a SingleCellExperiment object,
#' adding in feature information if provided and calculating QC metrics
#'
counts_to_sce <- function(input, separator = "\t", feature_data = "", cell_meta = list()) {
  data <- read.delim(input, sep = separator)
  sce <- SingleCellExperiment(assays = list(counts = data))
  db <- NULL
  if (feature_data && str_detect(feature_data, "\\.gtf")) {
    ensdb_file <- ensembledb::ensDbFromGtf(feature_data)
    db <- ensembldb::EnsDb(ensdb_file)
  } else if (feature_data && str_detect(feature_data, "\\.gff")) {
    ensdb_file <- ensembledb::ensDbFromGff(feature_data)
    db <- ensembldb::EnsDb(ensdb_file)
  } else if (feature_data && str_detect(feature_data, "\\.sqlite")) {
    db <- ensembldb::EnsDb(feature_data)
  }
  metric_subsets <- list()
  if (!is.null(db)) {
    sce <- add_feature_info(sce, db)
    metric_subsets$mito <- rowData(sce)$is_mito
  }
  for (c in cell_meta) {
    colData(sce)[[c]] <- c
  }
  scater::addPerCellQC(data, subsets = metric_subsets)
}

#' Wrapper function for calculating QC thresholds with MAD and determining outliers
#'  A cell is considered an outlier if it fails ANY of the metrics
#'
#' @param qc_spec a list of lists. keys are the name of the qc metric column to run the
#'    outlier detection on, and values are arguments passed to isOutlier for that metric
#' @param subfields character vector passed to sub.fields
#' @param batch_col column in the sce metadata containing batch information
#'
qc_mads <- function(sce, qc_spec = NULL, subfields = "subsets_mito_percent",
                    batch_col = NULL) {
  subfields <- keep(subfields, \(x) x %in% colnames(colData(sce)))
  batch <- NULL
  if (!is.null(batch_col)) batch <- sce[[batch_col]]
  reasons <- scuttle::perCellQCFilters(sce,
    sub.fields = subfields,
    batch = batch
  )
  if (!is.null(qc_spec)) {
    reasons_plus <- lapply(names(qc_spec), \(x) {
      result <- do.call(
        \(...) scuttle::isOutlier(sce[[x]], batch = sce[[batch_col]], ...),
        qc_spec[[x]]
      )
      to_df <- list()
      pars <- paste0(qc_spec[[x]], collapse = "_")
      colname <- paste0(x, "_", pars)
      to_df[[colname]] <- result
      df <- as.data.frame(to_df)
      df
    }) |>
      bind_cols()
    reasons <- cbind(reasons, reasons_plus)
    reasons$discard <- apply(reasons, 1, any)
  }
  thresholds <- lapply(colnames(reasons), \(x) {
    vec <- reasons[[x]]
    t <- attr(vec, "threshold")
    if (!is.null(t)) {
      tibble(metric = x, lower = t["lower"], higher = t["higher"])
    }
  }) |> bind_rows()
  qc_cols <- colnames(reasons) |> discard(\(x) x == "discard")
  list(df = reasons, thresholds = thresholds, qc_applied = qc_cols)
}

#' Aggregate count data into pseudobulks, returning
#' a DGEList
#'
#' @param cluster_col column of colData(sce) specifying cluster assignments used
#' for aggregation
#' @param method count aggregation method, one of "sum" or "mean"
#' @param char_agg function to aggregate colData() of type "character". By
#'  default takes unique values and collapses with comma delimiter
#' @param numeric_agg function to aggregate numeric colData(). Mean by default
#'
sce2pseudobulk <- function(sce, cluster_col = "cluster", pb_method = "sum",
                           aslot = "counts",
                           numeric_agg = \(x) mean(x, na.rm = TRUE),
                           char_agg = \(x) paste0(unique(sort(x)), collapse = ",")) {
  cluster <- colData(sce)[[cluster_col]]
  meta <- group_by(as_tibble(colData(sce)), !!as.symbol(cluster_col)) |>
    summarize(
      across(where(is.numeric), numeric_agg),
      across(where(is.character), char_agg)
    )

  if (pb_method == "sum") {
    cluster_mat <- model.matrix(~ 0 + cluster)
    aggregate <- assay(sce, aslot) %*% cluster_mat
  } else if (pb_method == "mean") {
    counts <- assay(sce, aslot)
    cells <- colnames(counts)
    aggregate <- do.call(cbind, lapply(unique(cluster), \(c) {
      current <- cells[cluster == c]
      rowMeans(counts[, current])
    }))
  }
  pb <- DGEList(aggregate, genes = rowData(sce), samples = meta)
  pb
}



#' Specify exactly for what reasons a cell is being discarded
#'
#' @description
#' Add a column to colData(sce) denoting the reason for which a cell is being
#' discarded in MADs-based filtering. By default, the full name of the reason is given
#' only if it is the sole reason - if a cell is discarded on 2+ reasons the number
#' is given instead
#'
identify_reasons <- function(sce, qc_cols, show_full = FALSE, reason_col = "discard_reason") {
  get_reason <- function(qc_row) {
    failed_qc <- purrr::keep(qc_row, \(x) x)
    if (length(failed_qc) == 0) {
      "kept"
    } else if (length(failed_qc) == 1) {
      names(failed_qc)
    } else if (show_full) {
      paste0(names(failed_qc), collapse = ";")
    } else {
      as.character(length(failed_qc))
    }
  }
  qc <- colData(sce)[, qc_cols]
  colData(sce)[[reason_col]] <- apply(qc, 1, \(x) get_reason(x))
  sce
}



#' Filter out cells with fixed thresholds
#'
#' @param sce_metrics a dataframe of sce metrics (row x metric) e.g. calculated
#' with perCellQCMetrics. Alternatively, an sce object where the metrics are in the
#' colData
#' @param threshold a list of metric->list(v=, d=),
#' where "metric" is a column in the metadata of the sce object
#' "d" (direction) is one of g|l|geq|leq or eq
#' "v" is the value
#' EX: if g, the cell value must be greater than the supplied value
qc_thresholds <- function(sce_metrics, thresholds) {
  threshold_tracker <- list()
  if (class(sce_metrics) == "SingleCellExperiment") {
    df <- colData(sce_metrics)
  }
  for (t in names(thresholds)) {
    if (t %in% colnames(df)) {
      val <- thresholds[[t]]$v
      direction <- thresholds[[t]]$d
      result <- switch(direction,
        g = {
          df[[t]] > val
        },
        l = {
          df[[t]] < val
        },
        geq = {
          df[[t]] >= val
        },
        leq = {
          df[[t]] <= val
        },
        eq = {
          df[[t]] == val
        }
      )
      threshold_tracker[[glue("{t}_{direction}_{val}")]] <- result
    }
  }
  df <- as.data.frame(threshold_tracker)
  df$discard <- apply(df, 1, any)
  df
}

#' Helper function for producing a diagnostic plot of cell loss
#'  Based on the plot described in Chapter 1.5 of OSCA Advanced
#'
#' @param gene_spec A list where names are cell types or conditions and values are the
#'    characteristic genes of that condition or type
#' @return A ggplot object plotting log-FC (discarded/kept) over log average count
diagnose_cell_loss <- function(sce, gene_spec, discard = NULL, identifier = "SYMBOL") {
  if (is.null(discard)) d <- sce$discard else d <- discard
  discarded <- calculateAverage(counts(sce)[, d])
  kept <- calculateAverage(counts(sce)[, !d])
  logged <- edgeR::cpm(cbind(discarded, kept), log = TRUE)
  id_vec <- rowData(sce)[[identifier]]
  log_fc <- logged[, 1] - logged[, 2] # Expression fold-change (discarded/kept)
  log_avg_counts <- rowMeans(logged)
  tb <- tibble(!!as.symbol(identifier) := id_vec,
    log_fc = log_fc,
    log_avg_counts = log_avg_counts
  )
  tb$Category <- map_chr(tb[[identifier]], \(x) {
    for (n in names(gene_spec)) {
      if (x %in% gene_spec[[n]]) {
        return(n)
      }
    }
    "Other"
  })
  ggplot(tb, aes(x = log_avg_counts, y = log_fc, color = Category)) +
    geom_point() +
    xlab("Log average count") +
    ylab("Log fold change (discarded/kept)")
}

plot_qc <- function(sce, y_axes, x_axis, titles = NULL, facet = NULL, discard_col = "discard") {
  plots <- list()
  if (!discard_col %in% colnames(colData(sce))) {
    discard_col <- NULL
  }
  if (is.null(titles)) {
    titles <- str_to_title(y_axes) |> str_replace_all("_", " ")
  }
  for (i in seq_along(y_axes)) {
    y <- y_axes[i]
    title <- titles[i]
    args <- list(x = x_axis, y = y, color_by = discard_col, show_violin = TRUE)
    if (!is.null(facet)) args$other_fields <- facet
    p <- do.call(\(...) plotColData(sce, ...), args) + ggtitle(title)
    if (!is.null(discard_col)) p <- p + guides(color = guide_legend(title = "Discard"))
    if (!is.null(facet)) p <- p + facet_wrap(as.formula(paste("~", facet)))
    if (i != 1 && !is.null(facet)) p <- p + theme(strip.text = element_blank(), strip.background = element_blank())
    if (i != length(y_axes)) {
      p <- p + theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
      )
    }
    plots[[i]] <- p
  }
  arranged <- ggarrange(plotlist = plots, common.legend = TRUE, ncol = 1)
  arranged
}

save_sce <- function(sce, filename) {
  if (str_detect(filename, "\\.hd5ad")) {
    zellkonverter::writeH5AD(sce, file = filename)
  } else {
    saveRDS(sce, filename)
  }
}

## * CLI entry point
scdbl_wrapper <- function(sce, args) {
  if (!is.null(args$scDblFinder_params)) {
    params <- rjson::fromJSON(file = args$scDblFinder_params)
    do.call(\(...) scDblFinder(sce, ...), args)
  } else {
    scDblFinder(sce)
  }
}


main <- function(args) {
  # TODO: <2025-01-23 Thu> must check that counts_to_sce works
  sce <- counts_to_sce(args$input)
  if (args$detect_doublets) {
    sce <- scdbl_wrapper(sce, args)
  }
  if (args$discard %in% c("both", "doublets")) {
    discarded <- sce[sce$scDblFinder.class != "singlet", ]
    sce <- sce[sce$scDblFinder.class == "singlet", ]
    outname <- paste0(args$output, "_doublets", ".", args$discard_suffix)
    save_sce(discarded, outname)
  }
  if (!is.null(args$simple_thresholds)) {
    thresh <- rjson::fromJSON(file = args$simple_thresholds)
    qc <- qc_thresholds(sce, thresh)
    colData(sce) <- cbind(colData(sce), qc)
  } else {
    qc_spec <- rjson::fromJSON(file = args$qc_spec)
    qc <- qc_mads(sce, qc_spec, subfields = args$subfields, batch_col = args$batch)
    colData(sce) <- cbind(colData(sce), qc$df)
    write_tsv(qc$thresholds, args$thresholds_output)
  }
  if (args$plot && !is.null(args$x_axis)) {
    y <- str_split_1(args$y_axes, ",")
    plot <- plot_qc(sce, y, args$x_axis, facet = args$facet)
    ggsave(filename = args$plot_name, plot = plot, dpi = 500, height = 20, width = 20)
  }
  if (!is.null(args$marker_genes)) {
    gs <- rjson::fromJSON(file = args$marker_genes)
    diagnose_cell_loss(sce, gene_spec = gs)
  }
  if (args$discard %in% c("both", "quality")) {
    discarded <- sce[sce$discard, ]
    sce <- sce[!sce$discard, ]
    outname <- paste0(args$output, "_low_quality", ".", args$discard_suffix)
    save_sce(discarded, outname)
  }
  save_sce(sce, args$output)
}


add_quality_args <- function(parser) {
  parser <- add_option(parser, c("-b", "--batch"),
    default = NULL,
    type = "character", help = "column containing batch data"
  )
  parser <- add_option(parser, "--subfields",
    default = "subsets_mito_percent",
    help = "Additional metric columns passed to the sub.fields argument of scuttle::perCellQCFilters"
  )
  parser <- add_option(parser, c("-m", "--n_mads"),
    default = 3, help = "Minimum number of MADs to consider a cell as an outlier"
  )
  parser <- add_option(parser, c("-q", "--qc_spec"),
    type = "character",
    help = "A json file defining additional qc parameters passed to scuttle::isOutlier.
Each key is the name of a column (metric) in the SingleCellExperiment object to operate on,
and the value is a list of arguments passed to isOutlier for that metric",
    default = NULL
  )
  parser <- add_option(parser, c("-t", "--thresholds_output"),
    type = "character",
    help = "Name of tsv file containing output thresholds from scuttle::perCellQCFilters",
    default = "thresholds.tsv"
  )
  parser <- add_option(parser, c("-s", "--simple_thresholds"),
    type = "character",
    help = "A json file defining the desired thresholds. If supplied, MAD thresholding is not carried out",
    default = NULL
  )
  parser <- add_option(parser, c("-p", "--plot"),
    type = "logical", default = TRUE,
    help = "Whether or not to produce diagnostic plots"
  )
  parser <- add_option(parser, c("-x", "--x_axis"),
    type = "character",
    help = "Metric to use as x axis",
    default = NULL
  )
  parser <- add_option(parser, c("-g", "--marker_genes"),
    type = "character",
    help = "path to a json file containing marker genes to use for diagnosing cell type loss",
    default = NULL
  )
  parser <- add_option(parser, c("-y", "--y_axes"),
    type = "character",
    help = "Comma-delimited list of SCE metric columns to plot on y axes",
    default = "sum,detected,subsets_mito_percent"
  )
  parser <- add_option(parser, c("-f", "--facet"), type = "character", help = "Facet of plot")
  parser <- add_option(parser, c("-n", "--plot_name"),
    type = "character",
    help = "Name of plot output file", default = "diagnostics.png"
  )
  parser <- add_option(parser, c("-l", "--loss_plot_name"),
    type = "character",
    help = "Name of loss plot output file", default = "loss_diagnostics.png"
  )
  parser
}

add_doublet_args <- function(parser) {
  parser <- add_option(parser, c("-d", "--detect_doublets"),
    type = "logical",
    help = "Whether or not to detect doublets"
  )
  parser <- add_option(parser, c("-r", "--scDblFinder_params"),
    type = "character",
    help = "json file of arguments to pass to scDblFinder::scDblFinder",
    default = NULL
  )
  parser
}

parse_args <- function() {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-d", "--discard"),
    type = "character",
    default = "none",
    action = "store_true",
    help = "Whether or not to remove cells based on the qc, or only flag them in metadata"
  )
  parser <- add_option(parser, "--discard_suffix",
    type = "character",
    help = "Suffix for file containing discarded cells", default = "h5ad"
  )
  parser <- add_option(parser, c("-i", "--input"), type = "character", help = "Input filename")
  parser <- add_option(parser, c("-o", "--output"), type = "character", help = "Output file name")
  parser <- add_quality_args(parser)
  parser <- add_doublet_args(parser)
  parse_args(parser)
}

if (sys.nframe() == 0) {
  args <- parse_args()
  main(args)
}
