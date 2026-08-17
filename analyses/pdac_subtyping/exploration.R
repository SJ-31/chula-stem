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

workdir <- here("analyses", "pdac_subtyping")
kras_status <- read_csv(here(workdir, "kras_status.csv")) |>
  mutate(status = factor(status))

ad <- import("anndata")
adata <- ad$read_h5ad(here(workdir, "cohort.h5ad"))
adata$obs <- merge(adata$obs, kras_status, on = "Case_ID")
adata <- adata[
  !is.na(adata$obs$status) & adata$obs$Project_ID == "CHULA_PHcase",
  !is.na(adata$var$GENENAME)
] # [2026-08-14 Fri]

filter_min_expr <- function(obj, min, n_samples = NULL) {
  if (!is.null(n_samples)) {
    keep <- rowSums(counts(obj) >= min) >= n_samples
  } else {
    keep <- rowSums(counts(obj)) >= min
  }
  obj[keep, ]
}

obj <- local({
  # do this for now
  obs <- adata$obs
  rownames(obs) <- obs$Case_ID
  x <- t(as.matrix(adata$X))
  storage.mode(x) <- "integer"
  rownames(x) <- rownames(adata$var)
  DESeqDataSetFromMatrix(
    countData = x,
    colData = obs,
    design = ~status
  )
})
obj <- filter_min_expr(obj, 10, 3)

dds <- DESeq(obj)
vsd <- varianceStabilizingTransformation(obj)
pca_plot <- plotPCA(vsd, "status")

counts_plot <- counts_violin(vsd)

rv <- rowVars(assay(vsd))
ntop <- 500
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca <- prcomp(t(assay(vsd)[select, ]))

clusters <- ConsensusClusterPlus(
  assay(vsd)[select, ]
)
# [2026-08-17 Mon] clustering of 2 had the best results
vsd$cluster <- clusters[[2]]$consensusClass[colnames(vsd)]

pca_clusters <- plotPCA(vsd, "status") +
  geom_text(aes(label = paste0(cluster, ": ", Case_ID)), color = "black")
# [2026-08-17 Mon] clusters coincide with splits from PCA
ggsave(
  here(workdir, "cases_pca.pdf"),
  plot = pca_clusters,
  width = 12,
  height = 9
)


greatest_vars <- pca$rotation |>
  as.data.frame() |>
  rownames_to_column() |>
  as_tibble() |>
  select(rowname, PC1) |>
  left_join(adata$var, by = join_by(x$rowname == y$GENEID)) |>
  arrange(PC1)

# Results for status 1 vs. 0
res <- results(dds)
res_shrunken <- lfcShrink(dds, res = res, coef = "status_1_vs_0")
res$log2FoldChange_shrink <- res_shrunken$log2FoldChange
res$lfcSE_shrink <- res_shrunken$log2FoldChange

res_tb <- as.data.frame(res) |>
  filter(padj <= 0.05) |>
  rownames_to_column() |>
  as_tibble() |>
  left_join(adata$var, by = join_by(x$rowname == y$GENEID)) |>
  arrange(desc(abs(log2FoldChange_shrink)))

res_tb |> write_csv(here(workdir, "de_results.csv"))

# Genes with large LFCs after shrinkage:
# Upregulated in KRAS positive: ROR2, TXNL4A
# Downregulated in KRAS positive: MYOM1, NGEF
## https://aacrjournals.org/cancerdiscovery/article/14/11/2162/749214/ROR2-Regulates-Cellular-Plasticity-in-Pancreatic
## https://www.nature.com/articles/srep12991

format_long <- function(obj, genes) {
  expr <- assay(obj)[genes, ] |>
    as.data.frame() |>
    rownames_to_column("gene_id") |>
    tidyr::pivot_longer(-gene_id, names_to = "sample")
  meta <- as.data.frame(colData(obj)) |> rownames_to_column("sample")
  dplyr::left_join(expr, meta, by = dplyr::join_by(sample))
}

long <- format_long(vsd, res_tb$rowname) |>
  left_join(adata$var, by = join_by(x$gene_id == y$GENEID)) |>
  left_join(res_tb, by = join_by(x$gene_id == y$rowname))

count_dist <- format_long(vsd) |>
  ggplot(aes(x = status, y = value)) +
  geom_violin()

norm_vals <- ggplot(
  long,
  aes(x = status, y = value, fill = log2FoldChange_shrink)
) +
  geom_boxplot() +
  facet_wrap(~GENENAME.x, scales = "free_y", ncol = 2) +
  guides(fill = guide_legend("Shrunken LFC")) +
  scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging")

ggsave(
  here(workdir, "de_norm_expr.pdf"),
  plot = norm_vals,
  width = 8,
  height = 20
)
