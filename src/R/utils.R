library(tidyverse)

read_with_filename <- function(x, col = "filename") {
  read_tsv(x) |> mutate(!!col := basename(x))
}

t2tb <- function(x, names = "rowname") {
  t(x) |>
    as.data.frame() |>
    rownames_to_column(var = names) |>
    as_tibble()
}

ensembl2entrez <- function(df) {
  mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  columns <- colnames(df)
  genes <- rownames(df)
  query <- biomaRt::getBM(
    attributes = c(
      "entrezgene_id", "ensembl_gene_id"
    ), filters = "ensembl_gene_id", values = genes,
    mart = mart
  )
  merged <- base::merge(df, query, by.x = 0, by.y = "ensembl_gene_id")
  merged <- merged[!is.na(merged$entrezgene_id), ] %>%
    distinct(entrezgene_id, .keep_all = TRUE)
  rownames(merged) <- merged$entrezgene_id
  merged[, columns]
}


basename_no_ext <- function(file) {
  bname <- basename(file)
  splits <- bname |> str_split_1("\\.")
  if (length(splits) > 1) {
    paste0(head(splits, n = -1), collapse = ".")
  } else {
    bname
  }
}

#' Flatten the character vector `char_vec` by `separator`, then
#'  get unique values
flatten_by <- function(char_vec, separator = ";", collapse = TRUE, unique = FALSE) {
  helper <- function(str) {
    if (is.na(str)) {
      return(NA)
    }
    lapply(str, \(x) {
      str_split_1(x, pattern = separator)
    }) |> unlist()
  }

  applied <- lapply(char_vec, helper) |> unlist()
  if (unique) {
    applied <- unique(applied)
  }
  if (collapse) {
    paste0(applied, collapse = separator)
  } else {
    list(applied)
  }
}

tb2map <- function(tb, keys, values, list = TRUE) {
  if (list) {
    as.list(tb[[values]]) |> `names<-`(tb[[keys]])
  } else {
    tb[[values]] |> `names<-`(tb[[keys]])
  }
}

into_char_list <- function(col, separator = ";") {
  lapply(col, \(x) {
    if (is.na(x)) {
      NA
    } else {
      str_split_1(x, ";")
    }
  })
}

modes <- function(x) {
  x <- discard(x, is.na)
  if (length(x) == 0) {
    return(NA)
  }
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

to <- function(obj, x, val) {
  obj[[x]] <- val
  obj
}


htest2tb <- function(test) {
  tibble(
    null = test$`null.value`,
    alternative = test$alternative,
    method = test$method,
    data = test$`data.name`,
    statistic = test$statistic,
    p_value = test$`p.value`
  )
}

get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


#' Retrieve RNA seq counts from separate count files and aggregate into
#' single genes X samples counts file
#'
#' @param metadata_tb a tb describing experiment metadata. Must have
#' at least a column `files` listing paths to all files and a column `cases` with
#' sample names
#' @param id_mapping a two-column tb or df where the first column is the old ids and
#'    the second column is the new
get_rnaseq_counts <- function(metadata_tb, id_mapping = NULL, sample_col = "cases",
                              file_col = "files") {
  sum_counts <- function(tb) {
    cols <- colnames(tb)
    gcol <- cols[1]
    scol <- cols[2]
    summed <- filter(tb, !is.na(!!as.symbol(gcol)) & !is.na(!!as.symbol(scol))) |>
      group_by(!!as.symbol(gcol)) |>
      summarise(!!as.symbol(scol) := sum(!!as.symbol(scol)))
    summed
  }

  if (!is.null(id_mapping)) {
    mapping_cols <- colnames(id_mapping)
    id_mapping <- dplyr::rename(id_mapping, gene_id = all_of(mapping_cols[1]))
    gcol <- mapping_cols[2]
  } else {
    gcol <- "gene_id"
  }
  counts <- apply(metadata_tb, 1, \(x) {
    col <- x[sample_col]
    tb <- read_tsv(x[file_col], col_names = c("gene_id", col))
    if (!is.null(id_mapping)) {
      joined <- inner_join(tb, id_mapping, by = join_by(gene_id))
      joined[, c(gcol, col)] |> sum_counts()
    } else {
      tb
    }
  }) |>
    reduce(\(x, y) left_join(x, y, by = join_by(!!as.symbol(gcol))))
  counts
}
