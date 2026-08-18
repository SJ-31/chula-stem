suppressMessages({
  options(
    reticulate.conda_binary = "/home/shannc/Bio_SDD/miniforge3/condabin/conda"
  )
  library(tidyverse)
  library(reticulate)
  library(pheatmap)
  library(DESeq2)
  library(paletteer)
  library(ConsensusClusterPlus)
  library(here)
  source(here("src", "R", "plotting.R"))
  use_condaenv("stem-base")
})

## * Utilities

format_long <- function(obj, genes) {
  expr <- assay(obj)[genes, ] |>
    as.data.frame() |>
    rownames_to_column("gene_id") |>
    tidyr::pivot_longer(-gene_id, names_to = "sample")
  meta <- as.data.frame(colData(obj)) |> rownames_to_column("sample")
  dplyr::left_join(expr, meta, by = dplyr::join_by(sample))
}

filter_min_expr <- function(obj, min, n_samples = NULL) {
  if (!is.null(n_samples)) {
    keep <- rowSums(counts(obj) >= min) >= n_samples
  } else {
    keep <- rowSums(counts(obj)) >= min
  }
  obj[keep, ]
}

prcomp_wrapper <- function(obj, ntop = 500) {
  rv <- rowVars(assay(obj))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  subset <- obj[select, ]
  pca <- prcomp(t(assay(subset)))
  greatest_vars <- merge(rowData(subset), pca$rotation, by = "row.names")
  list(pca = pca, vars = greatest_vars)
}


## * Design

## [2026-08-18 Tue] 'status' to check kras status,
## 'cluster' to check the weird clusterings
design_variable <- "status"
## design_variable <- "random"
## design_variable <- "cluster"
wanted_var_cols <- c("GENENAME", "GENEBIOTYPE")
workdir <- here("analyses", "pdac_subtyping")
rdir <- here(workdir, paste0("de_", design_variable))

## * Setup

dir.create(here(workdir, paste0("de_", design_variable)))
kras_status <- read_csv(here(workdir, "kras_status.csv")) |>
  mutate(status = factor(status))
cluster_status <- read_csv(here(workdir, "cluster_2026-08-18.csv")) |>
  mutate(status = factor(status)) |>
  rename(cluster = "status")

ad <- import("anndata")
adata <- ad$read_h5ad(here(workdir, "cohort.h5ad"))
adata$obs <- merge(adata$obs, kras_status, on = "Case_ID")
adata$obs <- merge(adata$obs, cluster_status, on = "Case_ID")
adata <- adata[
  !is.na(adata$obs$status) & adata$obs$Project_ID == "CHULA_PHcase",
  !is.na(adata$var$GENENAME)
] # [2026-08-14 Fri]

obj <- local({
  obs <- adata$obs
  rownames(obs) <- obs$Case_ID
  x <- t(as.matrix(adata$X))
  storage.mode(x) <- "integer"
  rownames(x) <- rownames(adata$var)
  DESeqDataSetFromMatrix(
    countData = x,
    colData = obs,
    design = as.formula(paste("~", design_variable)),
    rowData = adata$var,
    metadata = list(var_cols = colnames(adata$var))
  )
})
obj <- filter_min_expr(obj, 10, 3)

## * DE

dds <- DESeq(obj)
vsd <- varianceStabilizingTransformation(obj)
pca_plot <- plotPCA(vsd, "status")

counts_plot <- counts_violin(vsd, fill = "cluster") +
  theme(axis.text.x = element_text(angle = 90))

ggsave(
  here(workdir, "phcase_counts_dist_2026-08-18.png"),
  counts_plot,
  width = 10,
  height = 10
)
## clusters <- ConsensusClusterPlus(
##   assay(vsd)[select, ]
## )
# [2026-08-17 Mon] clustering of 2 had the best results
## vsd$cluster <- clusters[[2]]$consensusClass[colnames(vsd)]

pca_clusters <- plotPCA(vsd, "status", ntop = 500) +
  geom_text(aes(label = paste0(cluster, ": ", Case_ID)), color = "black")

# [2026-08-17 Mon] clusters coincide with splits from PCA
ggsave(
  here(workdir, "cases_pca.pdf"),
  plot = pca_clusters,
  width = 12,
  height = 9
)

pca_details <- prcomp_wrapper(vsd)
write_csv(pca_details$vars, here(workdir, "phcase_pca_greatest_vars.csv"))

included <- greatest_vars |>
  filter(GENEBIOTYPE != "protein_coding") |>
  pluck("rowname")

# Results for status 1 vs. 0
res <- results(dds)
res_shrunken <- lfcShrink(dds, res = res, coef = length(resultsNames(dds)))
res$log2FoldChange_shrink <- res_shrunken$log2FoldChange
res$lfcSE_shrink <- res_shrunken$log2FoldChange

res_tb <- as.data.frame(res) |>
  filter(padj <= 0.05) |>
  merge(rowData(dds)[, wanted_var_cols], by = "row.names") |>
  as.data.frame()

res_tb |>
  write_csv(here(workdir, paste0("de_", design_variable), "de_results.csv"))

v1 <- volcano_plot(res)
v2 <- volcano_plot(res, lfc_col = "log2FoldChange_shrink")
volcano <- v1 + v2 + plot_layout(guides = "collect")

# For design_variable == "status"
# Genes with large LFCs after shrinkage:
# Upregulated in KRAS positive: ROR2, TXNL4A
# Downregulated in KRAS positive: MYOM1, NGEF
## https://aacrjournals.org/cancerdiscovery/article/14/11/2162/749214/ROR2-Regulates-Cellular-Plasticity-in-Pancreatic
## https://www.nature.com/articles/srep12991

to_plot <- res_tb %>%
  as_tibble() %>%
  arrange(desc(abs(log2FoldChange_shrink))) %>%
  pluck("Row.names") %>%
  head(n = 25)


long <- format_long(vsd, to_plot) |>
  left_join(res_tb, by = join_by(x$gene_id == y$Row.names))

norm_vals <- ggplot(
  long,
  aes(x = status, y = value, fill = log2FoldChange_shrink)
) +
  geom_boxplot() +
  facet_wrap(~GENENAME, scales = "free_y", ncol = 2) +
  guides(fill = guide_legend("Shrunken LFC")) +
  scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging")

ggsave(
  here(workdir, paste0("de_", design_variable), "de_norm_expr.pdf"),
  plot = norm_vals,
  width = 8,
  height = 20
)

## * 2026-08-18 Weird clustering
## More evidence that the variation is due to technical differences
## - Large standard errors in DE analysis for genes with the highest
## fold changes
## - Direction of LFC is the same

res_tb |>
  group_by(GENEBIOTYPE) |>
  summarise(
    downreg = sum(log2FoldChange_shrink < 0),
    upreg = sum(log2FoldChange_shrink >= 0)
  ) |>
  write_csv(here(rdir, "de_direction_biotype.csv"))

# [2026-08-18 Tue] lncRNAs, TEC, miRNA are basically all downregulated
# in cluster 1 compared to 2

colData(vsd)$random <- sample(c("A", "B"), ncol(vsd), TRUE)
all_rna <- format_long(
  vsd
) |>
  left_join(res_tb, by = join_by(x$gene_id == y$Row.names)) |>
  filter(!is.na(GENEBIOTYPE)) |>
  mutate(biotype = GENEBIOTYPE) |>
  group_by(GENEBIOTYPE) |>
  mutate(count = n()) |>
  ungroup() |>
  mutate(GENEBIOTYPE = paste0(GENEBIOTYPE, "\nn = ", count, ""))
biotype_counts <- ggplot(
  all_rna,
  aes(
    x = !!as.symbol(design_variable),
    y = value,
    fill = !!as.symbol(design_variable)
  )
) +
  facet_wrap(~GENEBIOTYPE, scales = "free_y", ncol = 4) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = "none")
ggsave(
  here(rdir, "biotype_counts.pdf"),
  biotype_counts,
  width = 10,
  height = 15
)
group_by(all_rna, biotype, !!as.symbol(design_variable)) |>
  summarise(expr = round(mean(value), 3)) |>
  pivot_wider(names_from = !!as.symbol(design_variable), values_from = expr) |>
  write_csv(here(rdir, "biotype_means.csv"))
# TEC stands for To be Experimentally Confirmed

# [2026-08-18 Tue] Above observation can be confirmed with `biotype_counts`.
# Basically there is a systematic difference in expression for every
# biotype except protein coding genes
