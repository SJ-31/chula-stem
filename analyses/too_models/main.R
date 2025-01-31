library(here)
library(DGEobj.utils)
library(edgeR)
library(ensembldb)
library(glue)
library(tidyverse)
Sys.setenv(BIOMART_CACHE = here(".cache", "biomaRt"))
here::i_am("analyses/too_models/main.R")

U <- new.env()
source(here("src", "R", "utils.R"), local = U)

M <- list()
M$out <- here("analyses", "output", "too_models")
M$data <- here("analyses", "data")
M$remote <- here("analyses", "data_all")
M$id_mapping <- read_csv(here(M$data, "geneids_ensembl2entrez.csv"))
M$ensembl2entrez <- as.list(M$id_mapping$entrez) |> `names<-`(M$id_mapping$ensembl)
M$db <- EnsDb(here(M$data, "Homo_sapiens.GRCh38.113.sqlite"))
M$chula_raw_counts_file <- here(M$out, "chula_raw_counts.rds")
M$chula_tpm_file <- here(M$out, "chula_tpm.rds")
M$chula_meta_file <- here(M$out, "chula_metadata.tsv")

if (!file.exists(M$chula_raw_counts_file)) {
  tmp <- new.env()
  source(here("analyses", "too_models", "chula_data_prep.R"), local = tmp)
}
metrics <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")

chula_meta <- read_tsv(M$chula_meta_file)
chula_counts <- read_rds(M$chula_raw_counts_file) |>
  dplyr::filter(!gene_id %in% metrics) |>
  DGEList(samples = chula_meta)
chula_tpm <- read_rds(M$chula_tpm_file)
