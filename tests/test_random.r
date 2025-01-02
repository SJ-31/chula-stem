library(tidyverse)
library(GenomicRanges)
library(glue)
library(grid)
library(ggplot2)
cns_file <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.call.cns"
# Contains copy number calls
cnr_file <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.cnr"
# Contains the log2 segmentation data
cnn_file <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.targetcoverage.cnn"

facets_rds <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/5-patient_10-Facets/5-patient_10-Facets_hisens.rds"

format_chr_data <- function(tb, divide_by = 1000) {
  chr_levels <- c(as.character(1:22), "X", "Y") %>% paste0("chr", .)
  if (!str_detect(tb$chromosome[1], "^chr")) {
    tb$chromosome <- paste0("chr", tb$chromosome)
  }
  tb |> mutate(
    start = start / divide_by,
    end = end / divide_by,
    mid = (end - start) / 2 + start,
    chromosome = factor(tb$chromosome, levels = chr_levels)
  )
}

facets <- readRDS(facets_rds)
# Off-target bins typically going to be much wider than target bins
