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


prepare_counts <- function(dge, type) {
  if (type == "edgeR_cpm") {
    # Get cpm without normalization
    cur_counts <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  } else if (type == "edgeR_normalized_cpm") {
    # Get cpm with normalization, within the given tumor type i.e. batch
    cur_counts <- edgeR::cpm(edgeR::normLibSizes(dge),
      log = TRUE, prior.count = 1
    )
  } else if (type == "log2_fpkm") {
    gene_lengths <- AnnotationDbi::mapIds(M$db,
      keys = dge$genes[, 1],
      column = "SEQLENGTH", keytype = "GENEID"
    )
    cur_counts <- DGEobj.utils::convertCounts(dge$counts,
      unit = "fpkm", geneLength = gene_lengths, log = TRUE
    )
  } else {
    errorCondition("Type not recognized!")
  }
  cur_counts |>
    as.data.frame() |>
    mutate(gene_id = dge$genes[, 1]) |>
    relocate(gene_id, .before = everything()) |>
    as_tibble()
}

# <2025-01-30 Thu> Maybe add in normalized counts
quant_types <- c("edgeR_normalized_cpm", "kallisto_tpm", "edgeR_cpm", "log2_fpkm")
tumor_types <- unique(chula_meta$tumor_type)
naming_schemes <- list(
  entrez = \(x) {
    kept <- colnames(x)
    merged <- inner_join(x, M$id_mapping, by = join_by(x$gene_id == y$ensembl)) |>
      select(-gene_id) |>
      rename(gene_id = entrez) |>
      group_by(gene_id) |>
      mutate(across(where(is.numeric), mean)) |>
      ungroup() |>
      distinct(gene_id, .keep_all = TRUE)
    merged |> select(all_of(kept))
  },
  ensembl = NULL
)

# Produce a sample x genes csv, with a column for tumor_type
# TODO: make sure the models you're testing expects this data format
for (i in seq_along(tumor_types)) {
  t <- tumor_types[i]
  cur_type <- chula_meta |> filter(tumor_type == t)
  for (q in quant_types) {
    dir.create(here(M$out, q))
    if (q == "kallisto_tpm") {
      cur_counts <- chula_tpm |>
        select(gene_id, all_of(cur_type$cases)) |>
        mutate(across(where(is.numeric), log2))
    } else {
      cur <- chula_counts[, cur_type$cases]
      cur_counts <- prepare_counts(cur, q)
    }
    for (n in names(naming_schemes)) {
      outfile <- here(M$out, q, glue("chula-{q}-{t}-{n}.csv"))
      if (!is.null(naming_schemes[[n]])) {
        cur_counts <- naming_schemes[[n]](cur_counts)
      }
      transposed <- U$transpose(cur_counts, "gene_id") |> cbind(data.frame(tumor_type = t))
      write.csv(transposed, outfile)
    }
  }
}
