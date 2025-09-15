library(tidyverse)
library(glue)
library(paletteer)
library(SummarizedExperiment)

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
    rename = FALSE,
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

  if (rename) {
    tb[[response_col]] <- new_col
  } else {
    tb[[col_added]] <- new_col
  }
  tb
}

jaccard <- function(x, y) {
  length(intersect(x, y)) / length(union(x, y))
}

# Compute pairwise function `fn` between elements of `lst`
# if `symmetric`, then permutations of the same pair are computed
do_pairwise <- function(lst, fn, symmetric = FALSE) {
  tb <- tibble()
  names <- names(lst)
  for (j in seq_along(lst)) {
    for (k in seq_along(lst)) {
      if (k > j || symmetric) {
        row <- tibble_row(
          x = names[k],
          y = names[j],
          value = fn(lst[[k]], lst[[j]])
        )
        tb <- bind_rows(tb, row)
      }
    }
  }
  tb
}

## * Visualization

#' Produce visualizations and summary statistics for the output of
#' `response_group_ora`
#'
#' @param ora_results named list of tibble output from `response_group_ora`.
#'    names of the list correspond to the specific response `group_col` used in
#'    `response_group_ora`
#' @return A list containing...
#' 1. a list of venn diagrams, one for each response, showing the overlap of enriched
#'  terms between the labels of that group
#' 2. a single heatmap showing the Jaccard coefficient of the go terms between
#'  labels of each response. Specifically, between the groups in `response.label`
#'  e.g. `Paclitaxel.sensitive`
#' 3. If `individual`, and the ORA was done per-sample,
#'  a list of tibbles for each response denoting the percent enrichment
#'  of the GO term within the label
report_response_ora <- function(
    ora_results,
    sexp,
    individual = TRUE,
    palette = "ggthemes::Blue Light") {
  library(ggplot2)
  libray(ggVennDiagram)

  result <- list()

  helper <- function(group_col, ora_result) {
    result[[group_col]] <<- list()

    if (individual) {
      # Report the percentage of samples (per group label) that
      # ID was enriched in
      group_counts <- colData(exp) |>
        as_tibble() |>
        group_by(!!as.symbol(group_col)) |>
        summarise(group_size = n())

      tb <- ora_result |>
        group_by(ID, !!as.symbol(group_col)) |>
        summarise(n_samples = n()) |>
        left_join(group_counts, by = group_col) |>
        mutate(percentage = 100 * (n / group_size))

      result[[group_col]]$enrich_percent <<- distinct(
        ora_result,
        ID,
        Description
      ) |>
        inner_join(tb, by = join_by(ID))
    }

    grouped <- res |>
      group_by(!!as.symbol(group_col)) |>
      summarise(ID = list(ID))

    result[[group_col]]$venn <<- setNames(grouped$ID, grouped[[group_col]]) |>
      ggVennDiagram() +
      scale_fill_paletteer_c(palette)
  }

  all_responses <- lapply(names(ora_results), \(name) {
    helper(name, ora_results[[name]])
    ora_results[[name]] |>
      mutate(response_group = name) |>
      rename(label = name)
  }) |>
    bind_rows() |>
    mutate(new_lab = paste0(response_group, ".", label)) |>
    group_by(new_lab) |>
    summarise(ID = list(ID))

  response_list <- setNames(all_responses$ID, all_responses$new_lab)
  similarity <- do_pairwise(response_list, fn = jaccard, symmetric = FALSE)

  plot <- ggplot(similarity, aes(x = y, y = x, fill = value)) +
    geom_tile() +
    scale_x_discrete(position = "top") +
    xlab("ID") +
    ylab("ID") +
    scale_fill_manual(fill = "Jaccard") +
    scale_fill_paletteer_c(palette)
  result$jaccard <- plot

  result
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

## * Entry point

#' Main entry point for script
#' This function automatically writes results to `outdir` given on the cli
main <- function(exp, analysis, outdir, config) {
  if (analysis == "binary") {
    ## ** Binary testing
    responses <- config$response_cols
    for (resp in responses) {
      result <- binary_feature_analysis(exp, resp) |> arrange(p_adjust)
      write_tsv(result, glue("{outdir}/{resp}.tsv"))
    }
  } else if (analysis == "ora") {
    ## ** ORA
    responses <- colnames(colData(exp))

    setAnnotationHubOption("CACHE", config$ah_cache)
    ah <- AnnotationHub()
    orgdb <- query(ah, "org.Hs.eg.db")[[1]]

    background <- foo # TODO: get background
    result <- list()

    for (resp in responses) {
      tmp <- add_response_group(
        as_tibble(colData(exp)),
        response_col = resp,
        group_spec = config$response_groups,
        rename = TRUE
      )
      colData(exp) <- tmp
      result[[resp]] <- response_group_ora(
        exp,
        group_col = resp,
        orgdb = orgdb,
        individual = TRUE,
        background = background
      )
    }

    viz <- report_response_ora(
      result,
      exp,
      individual = TRUE,
      palette = config$palette_c
    )
  }
}

## * CLI

if (sys.nframe() == 0) {
  library(optparse)
  parser <- OptionParser()
  parser <- add_option(
    parser,
    c("-a", "--assay"),
    type = "character",
    help = "TSV file containing assay data. Columns are samples,
rows are features to associate with response"
  )
  parser <- add_option(
    parser,
    c("-r", "--analysis_routine"),
    type = "character",
    help = "Which analysis routine to run",
    default = NULL
  )
  parser <- add_option(
    parser,
    c("-c", "--config"),
    type = "character",
    help = "Config file",
    default = NULL
  )
  parser <- add_option(
    parser,
    c("-o", "--outdir"),
    type = "character",
    help = "Results output directory",
    default = NULL
  )
  args <- parse_args(parser)
  config <- jsonlite::read_json(args$config)
  exp <- make_exp(
    assay = dplyr::filter(read_tsv(args$assay), !is.na(symbol)),
    samples = read_tsv(config$response_data),
  )
  main(
    exp = exp,
    analysis = args$analysis_routine,
    config = config,
    outdir = args$outdir
  )
}
