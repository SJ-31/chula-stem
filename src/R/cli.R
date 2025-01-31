#!/usr/bin/env Rscript

library(tidyverse)
R_SRC <- Sys.getenv("R_SRC")
U <- new.env()
source(paste0(R_SRC, "/", "utils.R"), local = U)

if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-c", "--command"), type = "character", help = "Command to carry out")
  parser <- add_option(parser, c("-i", "--input"), type = "character", help = "Input filename")
  parser <- add_option(parser, c("-o", "--output"), type = "character", help = "Output file name")
  args <- parse_args(parser)
  if (args$command == "combine_counts") {
    spec <- read_csv(args$input)
    combined <- U$get_rnaseq_counts(spec, sample_col = "sample", file_col = "file")
    write_tsv(combined, args$output)
  }
}
