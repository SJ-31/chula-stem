R_SRC <- Sys.getenv("R_SRC")
source(paste0(R_SRC, "/", "utils.R"))
source(paste0(R_SRC, "/", "plotting.R"))
library(tidyverse)
library(here)
library(glue)
library(edgeR)
library(sva)

here::i_am("analyses/hcc/hccn_abundance.R")
data_dir <- here("analyses", "data_all")

format_counts <- function(tb) {
  filter(tb, !is.na(gene_name) & !is.na(count)) |>
    group_by(gene_name) |>
    summarise(count = sum(count))
}
outdir <- here("analyses", "output", "HCC_abundance")

save_fn <- function(plot, name) {
  ggsave(here(outdir, name), plot = plot, dpi = 500, width = 8, height = 8)
}
## * Data setup
counts_file <- here(outdir, "counts.rds")
if (!file.exists(counts_file)) {
  data_env <- new.env()
  source(here("analyses", "hcc", "hcc_data_setup.R"), local = data_env)
  rm(data_env)
} else {
  counts <- read_rds(counts_file)
}

M <- new.env()
M$tcgn_normals <- counts$samples |>
  filter(tissue == "N" & group == "TCGA") |>
  rownames()
M$wanted_sample <- "2-P17_tumor-STAR_Counts.tsv"
