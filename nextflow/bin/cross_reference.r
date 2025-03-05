#!/usr/bin/env Rscript
suppressMessages({
  library(tidyverse)
  library(GenomicRanges)
  R_SRC <- Sys.getenv("R_SRC")
  source(paste0(R_SRC, "/", "utils.R"))
})

##' Identify overlapping genomic ranges stored in tibbles x and y,
##'   for a specific chromosome and that are above a specific cutoff
##' @param x_id unique identifier column for ranges in x
##' @param y_id unique identifier column for ranges in y
##' @param chr column containing chromosome name, should be the same
##'   for both x and y
find_overlapping <- function(x, y, x_id, y_id,
                             x_start = "start", x_end = "end",
                             y_start = "start", y_end = "end",
                             chr_x = "chr", chr_y = "chr", percent_cutoff = 0.5) {
  x_g <- GRanges(x[[chr_x]], IRanges(x[[x_start]], x[[x_end]]), mcols = x[[x_id]])
  y_g <- GRanges(y[[chr_y]], IRanges(y[[y_start]], y[[y_end]]), mcols = y[[y_id]])
  pairs <- findOverlapPairs(x_g, y_g)
  overlaps <- pintersect(pairs)
  overlap_percent <- width(overlaps) / width(pairs@second)
  join_tb <- tibble(
    as.vector(pairs@first$mcols),
    as.vector(pairs@second$mcols),
    overlap_percent
  ) |>
    `colnames<-`(c(x_id, y_id, "p_overlap")) |>
    filter(p_overlap >= percent_cutoff)
  right_join(join_tb, x, by = {{ x_id }}) |> left_join(y, by = {{ y_id }})
}

overlapping_join <- function(
    x, y, x_on, y_on, x_start = "start", x_end = "end", y_start = "start",
    y_end = "end") {
  x_rename <- c(tempX = x_on, x_start = x_start, x_end = x_end)
  y_rename <- c(tempY = y_on, y_start = y_start, y_end = y_end)
  rename_back <- c(names(x_rename), names(y_rename)) |>
    `names<-`(c(x_rename, y_rename))
  x <- dplyr::rename(x, all_of(x_rename))
  y <- dplyr::rename(y, all_of(y_rename))
  inner_join(x, y, by = join_by(
    x$tempX == y$tempY,
    overlaps(x_start, x_end, y_start, y_end)
  )) |> rename(any_of(rename_back))
}

cross_reference_cnv <- function(input, reference, clingen, wc) {
  headers <- c(
    "VariantID", "Chromosome", "Type",
    "Classification", "Known or predicted dosage-sensitive genes",
    "All protein coding genes",
    "accession", "source", "ClinGen_report", "Start", "End", "Total", "score",
    "1A-B", "2A", "2B", "2C", "2D", "2E",
    "2F", "2G", "2H", "2I", "2J", "2K", "2L", "3",
    "4A", "4B", "4C", "4D", "4E", "4F-H", "4I", "4J", "4K",
    "4L", "4M", "4N", "4O", "5A", "5B", "5C", "5D", "5E", "5F", "5G", "5H", "p_overlap"
  )
  data <- read_tsv(input)
  if (nrow(data) == 0) {
    return(empty_tibble(headers))
  }
  data <- data |> filter(!is.na(Start) & !is.na(End))
  ref <- read_tsv(reference) |>
    filter(!is.na(start) & !is.na(end)) |>
    mutate(id = row_number())

  result <- find_overlapping(data, ref, "VariantID", "id", "Start", "End",
    chr_x = "Chromosome", percent_cutoff = 0.9
  ) |>
    filter(Type == type) |>
    select(all_of(c(colnames(data), wc)))

  if (!missing(clingen)) {
    dosage_data <- read_csv(clingen, comment = "+", skip = 4)
    # No need to perform filtering here because ClassifyCNV will have already determined
    # valid dosage-sensitive genes

    result <- result |>
      separate_longer_delim(`Known or predicted dosage-sensitive genes`, ",") |>
      left_join(
        dosage_data,
        by = join_by(x$`Known or predicted dosage-sensitive genes` == y$`GENE/REGION`), na_matches = "never"
      ) |>
      select(all_of(c(colnames(result), "ONLINE REPORT")))
  }
  result
}


cross_reference_msi <- function(input, reference, clingen, wc) {
  headers <- c(
    "VariantID", "chromosome", "gene_name",
    "left_flank_bases", "repeat_unit_bases", "right_flank_bases",
    "accession", "source", "ClinGen_report", "Start", "End",
    "gene_start", "gene_stop", "repeat_times", "difference",
    "P_value", "FDR", "rank", "p_overlap"
  )
  data <- read_tsv(input)
  if (nrow(data) == 0) {
    return(empty_tibble(headers))
  }
  data <- data |>
    filter(!is.na(start) & !is.na(end)) |>
    mutate(VariantID = row_number()) |>
    dplyr::rename(Start = "start", End = "stop")

  ref <- read_tsv(reference) |>
    filter(!is.na(start) & !is.na(end)) |>
    mutate(id = row_number(), chr = str_remove(chr, "^chr"))

  result <- find_overlapping(data, ref, "VariantID", "id",
    x_start = "Start", x_end = "End", chr_x = "chromosome",
    percent_cutoff = 0.01
  ) |> select(all_of(c(colnames(data), wc)))
  if (!missing(clingen)) {
    genes <- read_csv(clingen, comment = "+", skip = 4) |> filter(CLASSIFICATION == "Definitive")
    result <- left_join(result, genes,
      by = join_by(gene_name == `GENE SYMBOL`),
      na_matches = "never"
    ) |>
      select(all_of(c(colnames(result), "ONLINE REPORT")))
  }
  result
}

##' Cross-reference known CNVs or repetitive regions with those that were observed
##' @param type: either cnv or repeat|msi
##' @param input: if cnv data, output from ClassifyCNV.
##' If repeat data, a tsv file with at least the columns
##'     "chromosome, start, stop, gene_name"
##' @param clingen: for CNV data, a ClinGen dosage sensitivity data file
##'   for msi data, ClinGen gene-disease validity file
##' @return the original file with an added "Existing_variation" column containing
##' an accession or name of where that variant was found previously
cross_reference <- function(input, reference, type, clingen) {
  wanted_cols <- c("accession", "source", "p_overlap")
  if (str_to_lower(type) == "cnv") {
    result <- cross_reference_cnv(input, reference, clingen, wanted_cols)
  } else if (str_to_lower(type) == "msi") {
    result <- cross_reference_msi(input, reference, clingen, wanted_cols)
  }
  if (nrow(result) == 0) {
    return(result)
  }
  if (!missing(clingen)) {
    result <- dplyr::rename(result, ClinGen_report = "ONLINE REPORT")
  }
  result |>
    group_by(VariantID) |>
    filter(p_overlap == max(p_overlap) | is.na(p_overlap)) |>
    summarise(
      across(is.character, \(x) paste0(unique(discard(x, is.na)), collapse = ", ")),
      across(is.numeric, dplyr::first)
    ) |>
    mutate(across(is.character, \(x) na_if(x, "")))
}

if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-i", "--input"), type = "character", help = "Input filename")
  parser <- add_option(parser, c("-o", "--output"), type = "character", help = "Output file name")
  parser <- add_option(parser, c("-t", "--type"), type = "character", help = "Type of data (CNV|MSI)")
  parser <- add_option(parser, c("-r", "--reference"),
    type = "character",
    help = "Reference file containing regions of known variation"
  )
  parser <- add_option(parser, c("-c", "--clingen"),
    type = "character",
    help = "ClinGen dosage sensitivity file for CNV input, or gene-disease validity file for MSI input"
  )
  args <- parse_args(parser)
  result <- cross_reference(args$input, args$reference, args$type, args$clingen)
  write_tsv(result, args$output)
}
