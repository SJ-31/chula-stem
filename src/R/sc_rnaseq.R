R_SRC <- Sys.getenv("R_SRC")
source(paste0(R_SRC, "/", "utils.R"))
library(tidyverse)
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
#'
qc_wrapper <- function(sce, qc_spec = NULL, subfields = "subsets_mito_percent") {
  reasons <- scuttle::perCellQCFilters(sce, sub.fields = subfields)
  if (!is.null(qc_spec)) {
    reasons_plus <- lapply(names(qc_spec), \(x) {
      result <- do.call(\(...) scuttle::isOutlier(sce[[x]], ...), qc_spec[[x]])
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
