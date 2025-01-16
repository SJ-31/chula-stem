R_SRC <- Sys.getenv("R_SRC")
source(paste0(R_SRC, "/", "utils.R"))
source(paste0(R_SRC, "/", "plotting.R"))
library(tidyverse)
library(here)
library(glue)
library(edgeR)

here::i_am("analyses/hcc/hccn_abundance.R")

format_counts <- function(tb) {
  filter(tb, !is.na(gene_name) & !is.na(count)) |>
    group_by(gene_name) |>
    summarise(count = sum(count))
}

outdir <- here("analyses", "output", "HCC_abundance")
# <2025-01-15 Wed> Looking for upregulated expression of the MET gene (hepatocyte growth factor receptor)

## * Data setup
counts_file <- here(outdir, "counts.rds")
if (!file.exists(counts_file)) {
  data_env <- new.env()
  source(here("analyses", "hcc", "hcc_data_setup.R"), local = data_env)
  rm(data_env)
} else {
  counts <- read_rds(counts_file)
}

## * Perform DE

raw_plot <- pca_dgelist(counts, plot_aes = list(shape = "type", color = "group"))

normalized <- edgeR::normLibSizes(counts)

norm_plot <- pca_dgelist(normalized, plot_aes = list(shape = "type", color = "group"))
