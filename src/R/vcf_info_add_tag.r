#!/usr/bin/env Rscript

library(vcfR)
library(glue)

main <- function(name, description, number, type, default, input, output) {
  vcf <- read.vcfR(input, verbose = FALSE)
  info <- vcf@fix[, "INFO"]
  vcf@fix[, "INFO"] <- paste0(info, glue(";SOURCE={default}"))
  m <- glue('##INFO=<ID={name},Number={number},Type={type},Description="{description}">')
  vcf@meta <- append(vcf@meta, m)
  write.vcf(vcf, output)
  system2("gunzip", args = output)
}

if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-v", "--verbose"),
    action = "store_true",
    default = TRUE, help = "Print extra output [default]"
  )
  parser <- add_option(parser, c("-i", "--input"))
  parser <- add_option(parser, c("-o", "--output"),
    type = "character",
    help = "Output file name"
  )
  parser <- add_option(parser, c("-b", "--number"),
    type = "character", help = "Tag number"
  )
  parser <- add_option(parser, c("-t", "--type"), type = "character", help = "Type of value in tag")
  parser <- add_option(parser, c("-a", "--default"), type = "character", help = "Value of tag")
  parser <- add_option(parser, c("-d", "--description"), type = "character")
  parser <- add_option(parser, c("-n", "--name"), type = "character", help = "Name of tag")
  args <- parse_args(parser)
  output <- glue("{args$output}.gz")
  main(
    args$name, args$description, args$number,
    args$type, args$default, args$input, output
  )
}
