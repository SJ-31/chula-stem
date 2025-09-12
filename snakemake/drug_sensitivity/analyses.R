library(tidyverse)
library(glue)

## * Utils

#' Try to evaluate `expr`, returning a list of objects from
#' the attempt
#'
try_expr <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      err <<- e
      NULL
    }),
    warning = function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warning = warn, error = err)
}

#' Group samples according to specified response quantiles
#'
add_response_group <- function(
    tb,
    response_col,
    group_spec = list(
      sensitive = "<0.25",
      resistant = ">0.75",
      intermediate = NULL
    ),
    col_added = NULL) {
  if (is.null(col_added)) {
    col_added <- glue("{response_col}_group")
  }
  response <- tb[[response_col]]
  cdf <- ecdf(response)
  make_na <- is.na(response)
  to_check <- discard(group_spec, is.null)
  fill_group <- names(keep(group_spec, is.null))[1]

  checked <- lapply(names(to_check), \(g) {
    result <- list()
    spec <- group_spec[[g]]
    val <- cdf(response)
    vec <- eval(parse(text = glue("val{group_spec[[g]]}")))
    result[[g]] <- ifelse(vec, g, NA)
    as_tibble(result)
  }) |>
    bind_cols()
  new_col <- coalesce(!!!checked)
  new_col <- ifelse(is.na(new_col) & !make_na, fill_group, new_col)

  tb[[col_added]] <- new_col
  tb
}

## * Functions

#' Helper function for aligning assay data and sample metadata, creating
#' a summarized experiment object
#'
#' @param assay result of assay as a tibble. Columns are expected to be sample names, rows
#'  are features. One column should contain feature names
#' @param samples Sample metadata tibble i.e. colData
#' @param sample_col column in `samples` containing sample names
#'
make_exp <- function(
    assay,
    samples,
    join_on,
    how = "intersect",
    sample_col = "sample",
    assay_name = "x",
    assay_feature_col = "symbol") {
  row_data <- DataFrame(tibble(
    !!assay_feature_col := assay[[assay_feature_col]]
  ))
  assay_samples <- colnames(assay) |> discard(\(x) x == assay_feature_col)
  if (how == "intersect") {
    in_both <- intersect(samples[[sample_col]], colnames(assay))
    assay <- select(assay, all_of(c(assay_feature_col, in_both)))
    samples <- samples[samples[[sample_col]] %in% in_both, ]
  } else if (how == "keep_assay") {
    not_in <- assay_samples |> discard(\(x) x %in% samples[[sample_col]])
    samples <- samples[samples[[sample_col]] %in% assay_samples, ]
    extras <- tibble(!!sample_col := not_in)
    samples <- bind_rows(samples, extras)
  } else {
    stop(
      "A `how` method must be specified to combine the assay data and colData!"
    )
  }
  # Samples must be in the same order
  assay <- select(assay, all_of(c(assay_feature_col, samples[[sample_col]])))
  SummarizedExperiment(
    assays = `names<-`(
      list(column_to_rownames(
        assay,
        var = assay_feature_col
      )),
      assay_name
    ),
    colData = column_to_rownames(samples, var = sample_col),
    rowData = row_data
  )
}

#' For each variable (columns of `tb`), perform pairwise tests to
#'  see if there is a statistically significant difference in `response`
#'  between the subgroup of `tb` with and without the feature
#'
#' @param sexp SummarizedExperiment object containing assays and response data
#'  first assay assumed to be the one of interest
#' @param response response vector, should be aligned to samples in `tb`
#' @param test hypothesis testing function with two arguments. First argument is
#'      the response vector for the feature `present`
#' @param correction multiple-testing correction function, takes vector of p-values as argument
binary_feature_analysis <- function(
    sexp,
    response,
    test = \(x, y) wilcox.test(x, y, alternative = "greater"), # We are interested that the presence of a mutation
    correction = p.adjust) {
  response_vec <- colData(sexp)[[response]]
  features <- rownames(rowData(sexp))
  data <- assay(sexp)

  lapply(features, \(feat) {
    mask <- data[feat, ] > 0
    present <- response_vec[mask]
    absent <- response_vec[!mask]
    result <- list(feature = feat)
    test <- try_expr(test(present, absent))
    if (is.null(test$value)) {
      result$p_value <- NA
      result$success <- FALSE
      result$reason <- test$error$message
    } else {
      result$p_value <- test$value$p.value
      result$success <- TRUE
    }
    as_tibble(result)
  }) |>
    bind_rows() |>
    mutate(p_adjust = correction(p_value)) |>
    relocate(p_adjust, .after = p_value)
}

#' Identify pathways that are enriched in mutated genes for each group of
#' the factor `group_col`
#'
#' @param group_col Name of a factor column of rowData(sexp) to analyze
#' @param individual If True, enrich mutated genes for each sample of the subgroup
#' @param at_least if NOT individual, require this many hits of a gene in across the subgroup
#'    for the gene to be considered in the enrichment
response_group_ora <- function(
    sexp,
    group_col,
    orgdb,
    individual = TRUE,
    at_least = 1,
    background = NULL) {
  library(clusterProfiler)

  enrich_helper <- function(gene_list) {
    mapped_genes <- mapIds(
      orgdb,
      keys = gene_list,
      column = "ENTREZID",
      keytype = "SYMBOL"
    )
    as_tibble(enrichGO(gene = mapped_genes, orgdb, universe = background))
  }
  groupings <- colData(sexp)[[group_col]]
  lapply(unique(groupings), \(g) {
    mask <- groupings == g
    filtered <- sexp[, mask]
    if (individual) {
      sample_names <- colnames(filtered)
      result <- lapply(sample_names, \(x) {
        cur <- assay(filtered[, x])
        gene_list <- rownames(cur)[cur > 0]
        enrich_helper(gene_list) |> mutate(sample = x)
      }) |>
        bind_rows()
    } else {
      cur <- assay(filtered) |> rowSums()
      gene_list <- rownames(cur)[cur > 0]
      result <- enrich_helper(gene_list)
    }
    mutate(result, !!group_col := g)
  }) |>
    bind_rows()
}
