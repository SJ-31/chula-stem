R_SRC <- Sys.getenv("R_SRC")
source(paste0(R_SRC, "/", "utils.R"))
source(paste0(R_SRC, "/", "plotting.R"))
library(tidyverse)
library(here)
library(glue)
library(edgeR)
library(sva)

here::i_am("analyses/hcc/hccn_abundance.R")

format_counts <- function(tb) {
  filter(tb, !is.na(gene_name) & !is.na(count)) |>
    group_by(gene_name) |>
    summarise(count = sum(count))
}
outdir <- here("analyses", "output", "HCC_abundance")

save_fn <- function(plot, name) {
  ggsave(here(outdir, name), plot = plot, dpi = 500, width = 8, height = 8)
}

# <2025-01-15 Wed> Looking for upregulated expression of the MET gene (hepatocyte growth factor receptor)

## * Data setup
counts_file <- here(outdir, "counts.rds")
if (!file.exists(counts_file)) {
  data_env <- new.env()
  source(here("analyses", "hcc", "hcc_data_setup.R"), local = data_env)
  rm(data_env)
} else {
  counts <- read_rds(counts_file)
}

## * Batch effect visualization

to_keep <- filterByExpr(y = counts)
counts <- counts[to_keep, , keep.lib.sizes = FALSE]

to_keep <- filterByExpr(y = counts)
normalized <- edgeR::normLibSizes(counts)

norm_plot <- pca_dgelist(normalized, plot_aes = list(shape = "group", color = "type"))
raw_plot <- pca_dgelist(counts, plot_aes = list(shape = "group", color = "type"))
save_fn(norm_plot, "norm_pca.png")
save_fn(raw_plot, "raw_pca.png")

## ** Attempt correction with ComBat
full_model <- model.matrix(~ as.factor(tissue), data = counts$samples)

# This expects normalized and cleaned data
batch <- counts$samples$group
combat_corrected <- sva::ComBat(cpm(normalized, log = TRUE), batch = batch, mod = full_model)

combat_plot <- pca_dgelist(list(
  counts = combat_corrected,
  samples = normalized$samples
), plot_aes = list(shape = "group", color = "type"))
save_fn(combat_plot, "combat_pca.png")
# <2025-01-16 Thu> Looks good actually

# Combat version for raw count data
combat_seq_corrected <- sva::ComBat_seq(counts$counts,
  batch = batch,
  group = as.factor(counts$tissue),
  covar_mod = full_model
)
ccounts <- DGEList(counts = combat_seq_corrected, samples = counts$samples, genes = counts$genes) |>
  normLibSizes()
corrected_plot <- pca_dgelist(ccounts, plot_aes = list(shape = "group", color = "type"))
save_fn(corrected_plot, "combat_seq_pca.png")

## * DE
do_de <- function(dgelist) {
  data <- dgelist$samples
  design <- model.matrix(~ case + tissue, data = data)
  rownames(design) <- colnames(dgelist)
  disp <- estimateCommonDisp(dgelist, design)
  de <- exactTest(disp)
  result <- bind_cols(de$genes, de$table) |> as_tibble()
  list(de = result, disp = disp)
}

## chula_only
# <2025-01-16 Thu> Ordinarily you would include batch effects in the design here
# but in doing so, the treatment would only be compared within each batch and that won't
# work since the chula batch only has one treatment type (tumor)
# Might compare this to using the combat-corrected data
de_no_batch <- do_de(normalized)
de_corrected <- do_de(ccounts)

## * Interpret
de_no_batch$de |> filter(gene_name == "MET")
de_corrected$de |> filter(gene_name == "MET")
# <2025-01-16 Thu> MET is NOT differentially upregulated in the
# across tumors vs normal

# But what about for sample 17 specifically?
normalized_counts <- cpm(ccounts$counts) # CPM
tcga_normals <- local({
  names <- ccounts$samples |>
    filter(tissue == "N" & group == "TCGA") |>
    rownames()
  normalized_counts[, names] |> rowMeans()
})
wanted <- "2-P17_tumor-STAR_Counts.tsv"
merged <- bind_cols(tcga_normals, normalized_counts[, wanted], ccounts$genes) |>
  as_tibble() |>
  `colnames<-`(c("TCGA_normal", "P17", "gene_name"))
write_tsv(merged, here(outdir, "log2cpm_P17_tcga_normal.tsv"))
