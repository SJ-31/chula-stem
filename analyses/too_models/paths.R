library(here)
library(ensembldb)
library(tidyverse)
library(ensembldb)
library(DGEobj.utils)
library(edgeR)
library(glue)

out <- here("analyses", "output", "too_models")
data <- here("analyses", "data")
remote <- here("analyses", "data_all")
public <- here(remote, "public_data")
id_mapping <- read_csv(here(data, "geneids_ensembl2entrez.csv"))
ensembl2entrez <- as.list(id_mapping$entrez) |> `names<-`(id_mapping$ensembl)
db <- EnsDb(here(data, "Homo_sapiens.GRCh38.113.sqlite"))
chula_raw_counts_file <- here(out, "chula_raw_counts.rds")
chula_tpm_file <- here(out, "chula_tpm.rds")
chula_count_tpm_file <- here(out, "chula_tpm_scaled_count.rds")
chula_meta_file <- here(out, "chula_metadata.tsv")
tcga_data <- here(data, "tcga")

replace_ensembl_ids <- function(df, new_id_col) {
  kept <- colnames(df)
  merged <- inner_join(df, M$id_mapping, by = join_by(x$gene_id == y$ensembl)) |>
    dplyr::select(-gene_id) |>
    rename(c("gene_id" = new_id_col)) |>
    group_by(gene_id) |>
    mutate(across(where(is.numeric), mean)) |>
    ungroup() |>
    distinct(gene_id, .keep_all = TRUE)
  merged |> dplyr::select(all_of(kept))
}
