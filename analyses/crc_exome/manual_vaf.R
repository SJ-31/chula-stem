library(here)
library(tidyverse)
source(here("src", "R", "snakemake.R"))

tb <- read_tsv(here(
  "analyses",
  "output",
  "crc_exome",
  "2025-08-05",
  "vcfs",
  "combined.tsv"
))
config <- yaml::read_yaml(here(
  "analyses",
  "crc_exome",
  "vaf_config.yaml"
))

snakemake <- Snakemake(
  input = list(
    tmb = here(
      "analyses",
      "output",
      "crc_exome",
      "2025-08-06",
      "tumor_mutational_burden.tsv"
    ),
    sbs = here("analyses", "output", "crc_exome", "2025-08-06", "sbs.tsv"),
    combined_vep = here(
      "analyses",
      "output",
      "crc_exome",
      "2025-08-06",
      "vcfs",
      "combined.tsv"
    )
  ),
  params = list(wanted_genes = names(config$wanted_genes))
)
snakemake@config <- config
