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
  helper <- function(b) {
    splits <- b |> str_split_1("\\.")
    if (length(splits) > 1) {
      paste0(head(splits, n = -1), collapse = ".")
    } else {
      b
    }
  }
  map_chr(bname, helper)
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


read_existing <- function(filename, expr, read_fn = identity()) {
  if (file.exists(filename)) {
    read_fn(filename)
  } else {
    expr(filename)
  }
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
                              file_col = "files", gene_col = 1, count_col = 2,
                              read_fn = \(x) read_tsv(x, col_names = FALSE)) {
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
    tb <- read_fn(x[file_col]) |> dplyr::select(all_of(c(gene_col, count_col)))
    colnames(tb) <- c("gene_id", col)
    if (!is.null(id_mapping)) {
      joined <- inner_join(tb, id_mapping, by = join_by(gene_id))
      joined[, c(gcol, col)] |> sum_counts()
    } else {
      tb
    }
  }) |>
    reduce(\(x, y) full_join(x, y, by = join_by(!!as.symbol(gcol))))
  counts
}


#' Transpose a dataframe or tibble while explicitly specifying column names
#'
#' @param colnames either a vector of column names, or the index/name
#' of the new column names in the df
transpose <- function(df, colnames = 1) {
  if (length(colnames) != ncol(df) && is.numeric(colnames)) {
    tmp <- colnames
    colnames <- df[, colnames]
    df <- df[, -tmp]
  } else if (length(colnames) != ncol(df)) {
    tmp <- colnames
    colnames <- df[[colnames]]
    df[[tmp]] <- NULL
  }
  t(df) |> `colnames<-`(colnames)
}


x2counts <- function(sce) {
  assay(sce, "counts") <- assay(sce, "X")
  assay(sce, "X") <- NULL
  sce
}

#' Extracts a random subset of the dgelist, optionally by selecting from
#' a specific metadata column in the `samples` dataframe
#'
#' @param value take the random subset from samples with this value in `samples_col`
#' leaving the others unchanged
#' e.g. dgelist_random(X, 5, "tumor_type", value = "LIHC") takes
#'
dgelist_random <- function(dge, n, samples_col = NULL, value = NULL) {
  helper <- function(obj) {
    selection <- sample(colnames(obj), n)
    obj[, selection]
  }

  if (is.null(samples_col)) {
    helper(dge)
  } else if (!is.null(value)) {
    others <- dge[, dge$samples[[samples_col]]]
    chosen <- dge[, dge$samples[[samples_col]] == value]
    cbind(others, helper(chosen))
  } else {
    dges <- lapply(unique(dge$samples[[samples_col]]), \(x) {
      cur <- dge[, dge$samples[[samples_col]] == x]
      helper(cur)
    })
    purrr::reduce(dges, \(x, y) cbind(x, y))
  }
}


#' Parse a vcf location string of chr:start-end  in a tibble
#'
parse_loc <- function(tb, loc_col = "Loc", alt_col = "Alt", ref_col = "Ref") {
  separated <- separate_wider_delim(tb, all_of(loc_col),
    delim = ":", names = c("chromosome", "start")
  )
  if (!str_detect(separated$chromosome[1], "chr")) {
    separated$chromosome <- paste0("chr", separated$chromosome)
  }
  separated |>
    mutate(
      start = as.numeric(start),
      end = start + (nchar(!!as.symbol(alt_col)) - nchar(!!as.symbol(ref_col))),
      end = ifelse(end - start <= 0, start + 1, end)
    ) |>
    relocate(chromosome, start, end, .before = everything())
}

shift_stranded <- function(x, shift = 0L, ...) {
  GenomicRanges::shift(x, shift = shift * ifelse("-" == strand(x), -1, 1), ...)
}

into_granges <- function(vep_file,
                         allowed_consequences = c(
                           "missense_variant", "frameshift_variant",
                           "downstream_gene_variant", "upstream_gene_variant",
                           "stop_gained", "splice_region_variant", "inframe_deletion",
                           "splice_donor_5th_base_variant"
                         ),
                         wanted_cols = c(
                           "ref", "alt", "vaf", "alt_depth", "gene_biotype", "symbol",
                           "consequence", "existing_variation"
                         )) {
  library(GenomicRanges)
  if (is.character(vep_file)) {
    tb <- read_tsv(vep_file)
  } else {
    tb <- vep_file
  }

  tb <- tb |>
    separate_longer_delim(STRAND, ";") |>
    mutate(Consequence = lapply(Consequence, \(x) intersect(str_split_1(x, ";"), allowed_consequences))) |>
    dplyr::filter(map_lgl(Consequence, \(x) length(x) > 0)) |>
    rename_with(str_to_lower) |>
    dplyr::rename(gene_biotype = "biotype") |>
    mutate(strand = case_match(strand,
      "-1" ~ "-",
      "1" ~ "+",
      .default = "*"
    ))
  loc <- utils$parse_loc(tb, "loc", "alt", "ref")[, 1:3]

  gr <- GRanges(
    seqnames = loc$chromosome, ranges = IRanges(loc$start, loc$end),
    strand = Rle(tb$strand)
  )
  for (col in wanted_cols) {
    mcols(gr)[[col]] <- tb[[col]]
  }
  gr
}
