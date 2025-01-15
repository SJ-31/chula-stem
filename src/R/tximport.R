#!/usr/bin/env Rscript

R_SRC <- Sys.getenv("R_SRC")
source(paste0(R_SRC, "/", "utils.R"))
library(tximport)


get_tx2gene <- function(reference) {
  if (str_detect(reference, "\\.gtf|\\.gff")) {
    gtf <- rtracklayer::readGFF(reference) |> as_tibble()
  } else if (str_detect(reference, "\\.rds")) {
    gtf <- readRDS(reference)
  }
  tx2gene <- gtf |>
    dplyr::select(transcript_id, gene_id) |>
    dplyr::filter(!is.na(transcript_id))
  geneids2name <- gtf |>
    dplyr::select(gene_id, gene_name) |>
    unique()
  list(tx2gene = tx2gene, geneids2name = geneids2name)
}

import_one <- function(input, output, gtf, type, ignore_tx_version) {
  read_ref <- get_tx2gene(gtf)
  tx2gene <- read_ref$tx2gene
  geneids2name <- read_ref$geneids2name
  imported <- tximport(input,
    type = type, tx2gene = tx2gene,
    ignoreAfterBar = TRUE, ignoreTxVersion = ignore_tx_version
  ) |>
    as.data.frame() |>
    rownames_to_column(var = "gene_id") |>
    as_tibble() |>
    inner_join(geneids2name, by = join_by(gene_id))
  write_tsv(imported, output)
}

main <- function(input_directory, output_directory, gtf, prefix, suffix, type, filename_regex) {
  files <- list.files(input_directory, full.names = TRUE)
  ignore_tx_version <- FALSE
  read_ref <- get_tx2gene(gtf)
  tx2gene <- read_ref$tx2gene
  geneids2name <- read_ref$geneids2name
  lapply(files, \(x) {
    name <- basename_no_ext(x) %>% paste0(output_directory, "/", prefix, ., suffix)
    sample_tb <- read_tsv(x)
    imported <- tximport(x,
      type = type, tx2gene = tx2gene,
      ignoreAfterBar = TRUE, ignoreTxVersion = ignore_tx_version
    ) |>
      as.data.frame() |>
      rownames_to_column(var = "gene_id") |>
      as_tibble() |>
      inner_join(geneids2name, by = join_by(gene_id))
    write_tsv(imported, name)
  })
}

if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-i", "--input"), type = "character", help = "Input filename")
  parser <- add_option(parser, c("-o", "--output"), type = "character", help = "Output file name")
  parser <- add_option(parser, c("-t", "--type"),
    type = "character", help = "Type of input abundance file e.g. kallisto"
  )
  parser <- add_option(parser, c("-r", "--reference"),
    type = "character", help = "Reference file with transcript information (GTF/GFF or RDS)"
  )
  parser <- add_option(parser, c("-g", "--ignore_transcript_version"),
    type = "logical",
    action = "store_true",
    help = "Whether or not to ignore transcript version in reference file",
    default = FALSE
  )
  args <- parse_args(parser)
  import_one(
    args$input, args$output, args$reference, args$type,
    args$ignore_transcript_version
  )
}
