library(tidyverse)
library(ggplot2)
library(glue)
library(here)
library(reticulate)
library(cowplot)

here::i_am("./analyses/pdac/get_old_tsv.R")
utils <- new.env()
source(here("src", "R", "utils.R"), local = utils)

file <- here("analyses", "output", "pdac_vaf_merged_previous.rds")
data_path <- here("analyses", "data_all", "output", "PDAC", "OLD")
wanted_cols <- c(
  "CHROM", "POS", "REF", "ALT", "DP", "AD", "FILTER", "Effect", "dbSNP154_ID",
  "Gene_Name", "Feature_Type", "Feature_ID", "Transcript_BioType", "Putative_Impact"
)

if (!file.exists(file)) {
  files <- list.files(data_path, full.names = TRUE)
  tsvs <- lapply(files, \(x) {
    sample <- basename(x) |>
      str_remove("DAC") |>
      str_remove("\\.csv")
    tb <- read_csv(x) |>
      select(all_of(wanted_cols)) |>
      rename(SYMBOL = Gene_Name, Consequence = Effect, IMPACT = Putative_Impact) |>
      mutate(CHROM = str_remove(CHROM, "chr"), sample = sample)
    tb
  })
  merged <- bind_rows(tsvs)
  saveRDS(merged, file)
}
