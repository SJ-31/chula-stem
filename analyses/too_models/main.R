U <- new.env()
source(here::here("src", "R", "utils.R"), local = U)
P <- new.env()
source(here::here("src", "R", "plotting.R"), local = P)
M <- new.env()
source(here::here("analyses", "too_models", "paths.R"), local = M)

Sys.setenv(BIOMART_CACHE = here(".cache", "biomaRt"))

qn <- function(df) {
  cols <- colnames(df)
  rows <- rownames(df)
  data_sort <- apply(df, 2, sort)
  row_means <- rowMeans(data_sort)
  data_sort <- matrix(row_means,
    nrow = nrow(data_sort),
    ncol = ncol(data_sort),
    byrow = TRUE
  )
  index_rank <- apply(df, 2, order)
  normalized_data <- matrix(nrow = nrow(df), ncol = ncol(df))
  for (i in seq_len(ncol(df))) {
    normalized_data[, i] <- data_sort[index_rank[, i], i]
  }
  df <- as.data.frame(normalized_data)
  colnames(df) <- cols
  rownames(df) <- rows
  df
}

scale_counts <- function(dge) {
  dge$counts <- scale(dge$counts)
  dge
}

qn_counts <- function(dge) {
  dge$counts <- qn(dge$counts)
  dge
}


if (!file.exists(M$chula_raw_counts_file)) {
  tmp <- new.env()
  source(here("analyses", "too_models", "chula_data_prep.R"), local = tmp)
}
metrics <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")

chula_meta <- read_tsv(M$chula_meta_file)

chula_counts <- read_rds(M$chula_raw_counts_file) |>
  dplyr::filter(!gene_id %in% metrics) |>
  DGEList(samples = chula_meta, group = chula_meta$tumor_type)

tcga_counts <- local({
  dge <- readRDS(here(M$tcga_data, "tcga_all_tumors.rds"))
  dge <- dge[grepl("ENSG", dge$genes$gene_id), ]
  dge <- U$dgelist_random(dge, 35, samples_col = "tumor_type") # Balance number of tumor
  # types
  dge$genes$gene_id <- gsub("\\..*", "", dge$genes$gene_id)
  dge
})

tcga_meta <- as_tibble(tcga_counts$samples) |> mutate(cases = colnames(tcga_counts))

combined_counts <- local({
  common <- base::intersect(chula_counts$genes$gene_id, tcga_counts$genes$gene_id)
  c <- chula_counts[chula_counts$genes$gene_id %in% common, ]
  t <- tcga_counts[!duplicated(tcga_counts$genes), ]
  t <- t[t$genes$gene_id %in% common, ]
  t$samples$source <- "TCGA"
  c$samples$source <- "Chula"
  common_cols <- base::intersect(colnames(t$samples), colnames(c$samples))
  t$samples <- t$samples[, common_cols]
  c$samples <- c$samples[, common_cols]
  cbind(c, t)
})

plot <- P$pca_dgelist(scale_counts(combined_counts),
  plot_aes = list(color = "tumor_type", shape = "source"),
  to_cpm = FALSE
)
