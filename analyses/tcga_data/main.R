library(here)
library(edgeR)
library(glue)
library(tidyverse)
Sys.setenv(BIOMART_CACHE = here(".cache", "biomaRt"))

U <- new.env()
source(here("src", "R", "utils.R"), local = U)
P <- new.env()
source(here("src", "R", "plotting.R"), local = P)
remote <- here("analyses", "data_all")
outdir <- here("analyses", "data", "tcga")
public_data <- here(remote, "public_data")

get_tcga <- function(data_dir, tumor_type, strnd = "unstranded",
                     sequencer = "HiSeq2000") {
  manifest_file <- list.files(data_dir,
    pattern = "gdc_sample_sheet",
    full.names = TRUE
  )[1]
  sample_sheet <- read_tsv(manifest_file) |>
    rename_with(\(x) str_replace_all(x, " ", "_")) |>
    mutate(
      files = map2_chr(
        File_ID, File_Name,
        \(dir, file) paste(data_dir, dir, file, sep = "/")
      )
    )
  read_fn <- function(file) {
    read_tsv(file, skip = 1) |>
      select(gene_id, all_of(strnd)) |>
      filter(!is.na(gene_id))
  }

  counts <- U$get_rnaseq_counts(sample_sheet,
    read_fn = read_fn,
    sample_col = "Sample_ID"
  )

  DGEList(counts = counts, samples = sample_sheet)
}

rds_list <- list.files(outdir, pattern = "\\.rds", full.names = TRUE)
all_tcga <- lapply(rds_list, \(x) readRDS(x)) |> `names<-`(basename(rds_list))

combined <- local({
  dge <- purrr::reduce(all_tcga, \(x, y) {
    cbind(x, y)
  })
  dge$samples$tumor_type <- case_match(
    dge$samples$Project_ID,
    "TCGA-CHOL" ~ "CHOL",
    "TCGA-LIHC" ~ "LIHC",
    c("TCGA-COAD", "TCGA-READ") ~ "COADREAD"
  )
  dge$samples$group <- dge$samples$Project_ID
  dge$samples$Project_ID <- NULL
  dge
})


combined_tumors <- combined[, combined$samples$Sample_Type == "Primary Tumor"]
saveRDS(combined_tumors, here(outdir, "tcga_all_tumors.rds"))
