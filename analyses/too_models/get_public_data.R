M <- new.env()
source(here::here("analyses", "too_models", "paths.R"), local = M)
U <- new.env()
source(here("src", "R", "utils.R"), local = U)

## library("cBioPortalData")
## cBioPortalData::setCache(here(".cache", "cBioPortal"))
## cbio <- cBioPortal()

## studies <- getStudies(cbio, buildReport = TRUE)
## firehose <- filter(studies, grepl("Firehose", name))
## downloadStudy("lihc_tcga", use_cache = TRUE, ask = FALSE)

# HCC samples
others <- list(HCC = here(M$public, "GEO_GSE233421-hcc", "counts.csv"))

firehose_dirs <- list(
  BRCA = "TCGA_FIREHOSE_BREAST",
  LIHC = "TCGA_FIREHOSE_LIHC"
)

get_firehose <- function(f) {
  firehose_names <- names(firehose_dirs)
  firehose_meta <- list()
  firehose_all <- lapply(seq_along(firehose_dirs), \(x) {
    tb <- read_tsv(here(M$public, firehose_dirs[[x]], "data_mrna_seq_v2_rsem.txt"))
    cases <- tb |>
      dplyr::select(-Hugo_Symbol, -Entrez_Gene_Id) |>
      colnames()
    meta <- tibble(cases = cases, tumor_type = firehose_names[x])
    firehose_meta[[x]] <<- meta
    dplyr::select(tb, -Hugo_Symbol)
  })
  firehose_meta <- bind_rows(firehose_meta)
  firehose_all <- reduce(firehose_all, \(x, y) full_join(x, y, by = join_by(Entrez_Gene_Id))) |>
    distinct(Entrez_Gene_Id, .keep_all = TRUE) |>
    column_to_rownames(var = "Entrez_Gene_Id") |>
    mutate(across(where(is.numeric), \(x) replace_na(x, 0)))
  firehose_dge <- DGEList(firehose_all, samples = firehose_meta, genes = data.frame(gene_ids = rownames(firehose_all)))
  saveRDS(firehose_dge, f)
  firehose_dge
}

get_others <- function(f) {
  all <- lapply(names(others), \(x) {
    tb <- read_csv(others[[x]])
    cases <- tb |>
      select(where(is.numeric)) |>
      colnames()
    meta <- tibble(cases = cases, tumor_type = x)
    DGEList(tb, samples = meta, genes = data.frame(gene_ids = tb[, 1]))
  }) |>
    reduce(\(x, y) cbind(x, y))
  saveRDS(all, f)
  all
}

firehose <- U$read_existing(here(M$out, "firehose_counts.rds"), get_firehose, readRDS)
public <- U$read_existing(here(M$out, "public_counts.rds"), get_others, readRDS)
