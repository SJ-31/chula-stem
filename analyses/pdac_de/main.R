suppressMessages({
  library(tidyverse)
  library(here)
  library(DESeq2)
  library(patchwork)
  options(box.path = here("src"))
  suppressWarnings(box::use(R / utils[read_existing]))
})

workdir <- here("analyses", "pdac_de")
data <- here("analyses", "data_all")

tx2gene <- read_tsv(here("analyses", "data", "tx2gene.tsv"))

get_counts <- function(f) {
  library(reticulate)
  use_condaenv("stem-base")
  ad <- import("anndata")
  other_pdac <- ad$read_h5ad(here(
    data,
    "public_data",
    "h5ad",
    "in_house_organoids_salmon.h5ad"
  ))
  phcase <- read_tsv(here(
    data,
    "output",
    "PDAC-RNASEQ",
    "4-cohort-All_counts.tsv"
  )) |>
    filter(gene_id %in% other_pdac$var$gene_id)
  combined <- merge(
    as.data.frame(t(other_pdac$X)) |>
      `colnames<-`(other_pdac$obs$Case_ID) |>
      `rownames<-`(other_pdac$var$gene_id),
    column_to_rownames(phcase, var = "gene_id"),
    by = "row.names"
  ) |>
    merge(
      other_pdac$var,
      by.x = "Row.names",
      by.y = "gene_id"
    ) |>
    column_to_rownames("Row.names")
  write_csv(combined, f)
  combined
}


## * Data import

counts <- read_existing(here(workdir, "pdac_counts.csv"), get_counts, read_csv)
counts <- counts |>
  filter(!is.na(gene_name)) |>
  distinct(gene_name, .keep_all = TRUE) |>
  `rownames<-`(NULL) |>
  column_to_rownames("gene_name") |>
  rename_with(\(cols) {
    str_remove(cols, "_[12]_C") |> str_remove("_C")
  })
counts[is.na(counts)] <- 0
counts <- mutate(counts, across(everything(), as.integer))

pyrvinium <- local({
  tmp <- openxlsx::read.xlsx(
    here(
      workdir,
      "IC and Zscore PP update.xlsx"
    ),
    sheet = "IC50 80 90",
    skipEmptyRows = TRUE
  ) |>
    as_tibble() |>
    select(-c(X7, X8)) |>
    tail(n = -1) |>
    dplyr::rename(
      IC50 = "X2",
      IC80 = "X4",
      IC90 = "X6",
      X1 = "smallest.to.largest"
    )
  bind_rows(
    tibble(sample = tmp$X1, metric = "IC50", value = tmp$IC50),
    tibble(sample = tmp$X3, metric = "IC80", value = tmp$IC80),
    tibble(sample = tmp$X5, metric = "IC90", value = tmp$IC90)
  ) |>
    mutate(sample = str_remove(sample, " $")) |>
    pivot_wider(names_from = "metric", values_fn = max)
}) |>
  mutate(
    sample = str_trim(sample),
    sample = str_replace(sample, "PHcase", "PHcase_"),
    sample = str_replace(sample, " ", "_"),
    sample = str_remove(sample, "#3"),
    cohort = ifelse(str_starts(sample, "PHcase"), "PHcase", "PDAC")
  )


# Goal: determine relationship between gene expression and IC80 in PDAC

# TODO: read up on common ways of handling IC values

## * DE analysis

obs <- left_join(
  tibble(sample = colnames(counts)),
  pyrvinium,
  by = join_by(sample)
) |>
  filter(!is.na(IC80)) |>
  mutate(across(starts_with("IC"), as.double))

# [2026-07-15 Wed] NOTE: taking logarithm of IC80 doesn't change anything
# TODO: experiment with

dds <- DESeqDataSetFromMatrix(
  countData = counts[, obs$sample],
  colData = column_to_rownames(obs, "sample"),
  design = ~ 0 + cohort + IC80
)
dds <- dds[rowSums(counts >= 10) >= 5, ]
dds <- DESeq(dds)

# Important to note that non-PHcase samples were prepared with
# salmon counts

ic80_res <- results(dds, name = "IC80")

## * Visualizations
# normalized data
vsd <- vst(dds, blind = FALSE)
cohort_plot <- plotPCA(vsd, intgroup = c("cohort"))
# Batch separation is clear

## plotPCA(vsd, intgroup = c("IC80"))

bplot <- obs |>
  mutate(log2_IC80 = log2(IC80)) |>
  select(sample, IC80, log2_IC80) |>
  pivot_longer(-sample) |>
  ggplot(aes(y = sample, x = value)) +
  geom_bar(stat = "identity") +
  facet_wrap(~name, scales = "free_x")
