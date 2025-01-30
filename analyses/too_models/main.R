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
M$ensembl2entrez <- read_csv(here(M$data, "geneids_ensembl2entrez.csv"))
M$db <- EnsDb(here(M$data, "Homo_sapiens.GRCh38.113.sqlite"))

## mapIds(db, keys = )

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
  filter(!gene_id %in% metrics) |>
  DGEList()
chula_tpm <- read_rds(M$chula_tpm_file)

dge <- DGEList(chula_counts)

q()

# <2025-01-30 Thu> Maybe add in normalized counts
quant_types <- c("edgeR_normalized_cpm", "kallisto_tpm", "edgeR_cpm", "fpkm")
tumor_types <- unique(chula_meta$tumor_type)
for (i in seq_along(tumor_types)) {
  t <- tumor_types[i]
  cur_type <- chula_meta |> filter(tumor_type == t)
  for (q in quant_types) {
    outfile <- here(M$out, q, glue("chula-{q}-{t}.csv"))
    if (q == "kallisto_tpm") {
      cur_counts <- chula_tpm |>
        select(gene_id, all_of(cur_type$cases)) |>
        mutate(tumor_type = t) |>
        mutate(across(where(is.numeric), log2))
    } else {
      cur <- chula_counts[, cur_type$cases]
      if (q == "edgeR_cpm") {
        # Get cpm without normalization
        cur_counts <- edgeR::cpm(cur, log = TRUE, prior.count = 1) |>
          as.data.frame() |>
          cbind(data.frame(tumor_type = t))
      } else if (q == "edgeR_normalized_cpm") {
        # Get cpm with normalization, within the given tumor type i.e. batch
        cur_counts <- edgeR::cpm(edgeR::normLibSizes(cur),
          log = TRUE, prior.count = 1
        ) |>
          as.data.frame() |>
          cbind(data.frame(tumor_type = t))
      } else if (q == "log2_fpkm") {
        gene_lengths <- AnnotationDbi::mapIds(M$db,
          keys = cur$genes[, 1],
          column = "SEQLENGTH", keytype = "GENEID"
        )
        cur_counts <- DGEobj.utils::convertCounts(cur,
          unit = "fpkm", geneLength = gene_lengths, log = TRUE
        ) |>
          as.data.frame() |>
          cbind(data.frame(tumor_type = t))
      }
    }
    write_tsv(cur_counts, outfile)
  }
}
