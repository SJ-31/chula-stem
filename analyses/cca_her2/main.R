suppressMessages({
  library(tidyverse)
  library(DESeq2)
  library(here)
  library(paletteer)
  library(patchwork)
  options(box.path = here("src"))
  suppressWarnings(box::use(R / utils[alr_normalize, swap_ids]))
})

# Since tumor type will be confounded with batch, you can't do a proper
# DE analysis. Instead just calculate ratios of TPM expression of the
# ERBB2 gene to housekeeping genes and compare across samples

id_mapping <- read_csv(here("analyses", "data", "ensembl_gene_data.csv")) |>
  select(ensembl_gene_id, hgnc_symbol)

workdir <- here("analyses", "cca_herr2")
remote <- here("analyses", "data_all")
brca_data <- here("analyses", "data", "brca_metadata")

# BRCA data sources
# - Early on-treatment transcriptional profiling as a tool for improving pathological response prediction in HER2-positive inflammatory breast cancer https://doi.org/10.1177/17588359221113269
# - Integrated DNA and RNA Sequencing Reveals Drivers of Endocrine Resistance in Estrogen Receptor–Positive Breast Cancer Open Access https://doi.org/10.1158/1078-0432.CCR-21-3189

## * Prepare data

cca_raw <- here(remote, "output/CCA/RNASEQ/4-cohort-All_counts.tsv") |>
  read_tsv() |>
  swap_ids(id_mapping, "gene_id", "hgnc_symbol", "ensembl_gene_id")

metadata <- tibble(
  sample = colnames(cca_raw),
  tumor_type = "CCA",
  source = "Chula"
)

ov_samples <- read_csv(here(
  brca_data,
  "overmoyer2021",
  "metadata_corrected.csv"
)) |>
  filter(treatment == "pre") |>
  pluck("samplename")

metadata <- bind_rows(
  metadata,
  tibble(sample = ov_samples, tumor_type = "BRCA", source = "overmoyer2021")
)

overmoyer2021 <- read_csv(here(
  brca_data,
  "overmoyer2021",
  "counts.raw.csv.gz"
)) |>
  swap_ids(id_mapping, "ensembl_gene_id", "hgnc_symbol")
overmoyer2021 <- overmoyer2021[, ov_samples]

p_samples <- read_tsv(here(brca_data, "perez2022", "E-MTAB-9917.sdrf.txt")) |>
  inner_join(
    read_tsv(here(
      brca_data,
      "perez2022",
      "patient_characteristics.tsv"
    )),
    by = join_by(x$`Characteristics[individual]` == y$`Patient ID`)
  ) |>
  filter(`HER2 Status` == "pos") |>
  group_by(`Characteristics[individual]`) |>
  summarise(source = dplyr::first(sort(`Source Name`))) |>
  pluck("source")

metadata <- bind_rows(
  metadata,
  tibble(sample = p_samples, tumor_type = "BRCA", source = "perez2022")
)

perez2022 <- read_csv(here(brca_data, "perez2022", "counts.csv")) |>
  column_to_rownames(var = "...1")
perez2022 <- perez2022[, p_samples]

obj <- DESeqDataSetFromMatrix(
  purrr::reduce(
    list(perez2022, cca_raw, overmoyer2021),
    \(x, y) column_to_rownames(merge(x, y, by = "row.names"), "Row.names")
  )[, metadata$sample] |>
    mutate(across(everything(), as.integer)),
  colData = column_to_rownames(metadata, var = "sample"),
  design = ~ 0 + tumor_type
)
obj <- estimateSizeFactors(obj)


## * Viz

h2 <- c("FOXH1", "CD2BP2", "KRT8P33", "RHOQP1", "GPI", "AUH", "ATG3", "ARF1")


her2 <- alr_normalize(obj, "ERBB2", h2) |>
  merge(colData(obj), by.x = "sample", by.y = "row.names") |>
  as.data.frame() |>
  inner_join(
    tibble(
      sample = colnames(obj),
      counts = unlist(counts(obj, normalize = TRUE)["ERBB2", ])
    ),
    by = join_by(sample)
  )


ratio_plot <- ggplot(her2, aes(x = reference, y = n_expr, fill = tumor_type)) +
  geom_boxplot() +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)
  ) +
  ylab("log(ratio of housekeeping gene expression)") +
  xlab("Housekeeping gene") +
  guides(
    fill = guide_legend("Tumor type")
  ) +
  scale_fill_paletteer_d("ggthemes::wsj_rgby")

n_counts <- ggplot(
  filter(her2, reference == head(her2$reference, n = 1)),
  aes(x = tumor_type, y = counts, fill = tumor_type)
) +
  geom_boxplot() +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)
  ) +
  ylab("Normalized expression") +
  guides(
    fill = guide_legend("Tumor type")
  ) +
  scale_fill_paletteer_d("ggthemes::wsj_rgby") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ggtitle("HER2 Expression Comparison")

plot <- n_counts +
  ratio_plot +
  plot_layout(guides = "collect", widths = c(1, 5))

ggsave(here(workdir, "plot.png"), width = 16, height = 10, dpi = 500)
