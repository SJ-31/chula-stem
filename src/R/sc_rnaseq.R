R_SRC <- Sys.getenv("R_SRC")
source(paste0(R_SRC, "/", "utils.R"))
library(edgeR)
library(tidyverse)
library(glue)
library(zellkonverter)

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
    anno_cols <- c("SYMBOL", "GENEID", "GENEBIOTYPE", "SEQNAME")
    anno <- AnnotationDbi::select(db,
      keys = rownames(data), keytype = "GENEID",
      columns = anno_cols
    )
    rowData(sce) <- left_join(sce.frame(GENEID = rownames(sce)), anno,
      by = join_by(GENEID)
    )
    is_mito <- (rowData(data)$SEQNAME == "MT") %>% replace_na(FALSE)
    metric_subsets$mito <- is_mito
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
qc_wrapper <- function(sce, qc_spec = NULL, subfields = "subsets_mito_percent",
                       batch_col = NULL) {
  reasons <- scuttle::perCellQCFilters(sce, sub.fields = subfields, batch = batch_col)
  if (!is.null(qc_spec)) {
    reasons_plus <- lapply(names(qc_spec), \(x) {
      result <- do.call(
        \(...) scuttle::isOutlier(sce[[x]], batch = batch_col, ...),
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
  list(df = reasons, thresholds = thresholds)
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


main <- function() {

}
