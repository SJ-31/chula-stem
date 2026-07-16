suppressMessages({
  library(clusterProfiler)
  library(DESeq2)
  library(tidyverse)
  library(here)
  library(glue)
  library(paletteer)
  library(patchwork)
  library(org.Hs.eg.db)
  options(box.path = here("src"))
  suppressWarnings(box::use(R / utils[read_existing], dplyr[select, filter]))
  set.seed(240)
})


workdir <- here("analyses", "pdac_de")
data <- here("analyses", "data_all")

hallmark <- read_csv(
  here("analyses", "data", "hallmark2symbol.csv"),
  col_names = c("term", "symbol")
)
reactome <- read_csv(
  here("analyses", "data", "reactome2symbol.csv"),
  col_names = c("term", "symbol")
)

## tx2gene <- read_tsv(here("analyses", "data", "tx2gene.tsv"))

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
# TODO: check if you should scale and center

## * Analysis

## ** DE

obs <- left_join(
  tibble(sample = colnames(counts)),
  pyrvinium,
  by = join_by(sample)
) |>
  filter(!is.na(IC80)) |>
  mutate(across(starts_with("IC"), as.double)) |>
  filter(IC80 > 0)


dds <- DESeqDataSetFromMatrix(
  countData = counts[, obs$sample],
  colData = column_to_rownames(obs, "sample"),
  design = ~ 0 + cohort + scale(log(IC80))
)
dds <- dds[rowSums(counts >= 10) >= 5, ]
dds <- DESeq(dds)

# Important to note that non-PHcase samples were prepared with
# salmon counts

# Fold-change of the covariate is the per-unit increase in the covariate
ic80_res <- results(dds)
ic80_res <- ic80_res[ic80_res$padj <= 0.05, ] |> as.data.frame()

ic80_res |>
  rownames_to_column(var = "symbol") |>
  write_csv(here(workdir, "de_genes.csv"))

## ** Enrichment

enrich_list <- ic80_res$log2FoldChange |>
  `names<-`(rownames(ic80_res)) |>
  sort(decreasing = TRUE)

## gse_res <- gseGO(
##   gene = enrich_list,
##   OrgDb = org.Hs.eg.db,
##   keyType = "SYMBOL",
##   ont = "ALL"
## )
## min(gse_res$p.adjust, na.rm = TRUE) # 0.1725709
## [2026-07-16 Thu] No significant terms identified with gsea

res_list <- list(
  downreg = ic80_res[ic80_res$log2FoldChange < 0, ],
  upreg = ic80_res[ic80_res$log2FoldChange > 0, ]
)

do_enrichment <- function(f) {
  tb <- lapply(names(res_list), \(g) {
    genes <- rownames(res_list[[g]])
    ora_go <- enrichGO(
      gene = genes,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "ALL",
      universe = rownames(dds)
    ) |>
      as_tibble() |>
      mutate(group = paste0("GO:", ONTOLOGY)) |>
      select(-ONTOLOGY)
    ora_hallmark <- enricher(
      genes,
      TERM2GENE = hallmark,
      universe = rownames(dds)
    ) |>
      as_tibble() |>
      mutate(group = "hallmark")
    ora_reactome <- enricher(
      genes,
      TERM2GENE = reactome,
      universe = rownames(dds)
    ) |>
      as_tibble() |>
      mutate(group = "Reactome")
    bind_rows(ora_hallmark, ora_reactome, ora_go) |> mutate(Direction = g)
  }) |>
    bind_rows()
  write_csv(tb, f)
  tb
}

enrich_res <- read_existing(
  here(workdir, "enrichment.csv"),
  do_enrichment,
  read_csv
)

## * Visualizations
# normalized data
vsd <- vst(dds, blind = FALSE)
vsd$log_IC80 <- log(vsd$IC80)

cohort_plot <- plotPCA(vsd[rownames(ic80_res), ], intgroup = c("cohort")) +
  guides(color = guide_legend("cohort"))
ic80_plot <- plotPCA(vsd[rownames(ic80_res), ], intgroup = c("log_IC80")) +
  scale_color_paletteer_c("ggthemes::Green-Gold") +
  guides(color = guide_legend("log(IC80)")) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())

pca_plot <- cohort_plot +
  ic80_plot +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "PCA of PDAC samples",
    subtitle = "After subsetting by DE genes"
  )
ggsave(
  filename = here(workdir, "pca_plot.png"),
  plot = pca_plot,
  height = 10,
  width = 10,
  dpi = 500
)

ic80_heatmap <- as.data.frame(assay(vsd)[rownames(ic80_res), ]) |>
  `colnames<-`(colnames(vsd)) |>
  rownames_to_column(var = "gene") |>
  as_tibble() |>
  pivot_longer(-gene, names_to = "sample", values_to = "expr") |>
  inner_join(rownames_to_column(ic80_res, "gene"), by = join_by(gene)) |>
  inner_join(
    rownames_to_column(as.data.frame(colData(vsd)), "sample"),
    by = join_by(sample)
  ) |>
  ggplot(
    aes(
      x = cut(log(IC80), breaks = 10),
      y = factor(gene, levels = names(enrich_list)),
      fill = expr
    )
  ) +
  geom_tile(width = 0.98) +
  theme_minimal() +
  scale_fill_paletteer_c("grDevices::Sunset") +
  guides(fill = guide_legend("Normalized expression")) +
  xlab("Binned log(IC80)") +
  ylab("Gene")
ggsave(
  plot = ic80_heatmap,
  filename = here(workdir, "ic80_heatmap.pdf"),
  height = 10,
  width = 15
)

lapply(rownames(ic80_res), \(gene) {
  lfc <- ic80_res[gene, ]["log2FoldChange"] |> round(3)
  p_val <- ic80_res[gene, ]["padj"] |> round(4)
  plot <- tibble(
    expr = assay(vsd)[gene, ],
    IC80 = colData(vsd)$IC80,
    cohort = colData(vsd)$cohort
  ) |>
    ggplot(aes(x = log(IC80), y = expr, color = cohort)) +
    geom_point() +
    ylab("Normalized expression") +
    labs(
      title = gene,
      subtitle = glue("log2FoldChange: {lfc}, p-value: {p_val}")
    )
  ggsave(filename = here(workdir, glue("gene_plots/{gene}.pdf")), plot = plot)
})

## * GATA6

# TODO: Gata6 level, normalized against housekeeping genes. In pdac
# Just make it a bar plot and facet by the housekeeping gene that was normalized
# Use ALR on the raw counts
housekeeping <- c("MLF2", "SF3B1", "GNB1", "CTBP1", "MYH9")
# Took the top 5 genes in pancreas from HRT atlas, sorting by lowest STD

normalize_housekeeping <- function(obj, target, housekeeping) {
  lapply(housekeeping, \(hk) {
    normalized <- log(assay(obj)[target, ], assay(obj)[hk, ])
    tibble(
      sample = names(normalized),
      n_expr = normalized,
      housekeeping = hk
    ) |>
      mutate(rank = rank(-n_expr))
  }) |>
    bind_rows()
}

h2 <- c("FOXH1", "TMEM269", "CD2BP2", "KRT8P33", "RHOQP1", "FREY1")

gata6_hk_plot <- normalize_housekeeping(
  dds[, str_detect(colnames(dds), "PHcase")],
  housekeeping = housekeeping,
  target = "GATA6"
) |>
  ggplot(aes(x = sample, y = housekeeping, fill = n_expr)) +
  geom_tile() +
  geom_text(aes(label = rank), color = "white", size = 6) +
  xlab("Sample") +
  ylab("Housekeeping gene") +
  guides(fill = guide_legend("Normalized expression\n(additive log ratio)")) +
  labs(
    title = "GATA6 Expression"
  )
gata6_hk_plot

ggsave(
  here(workdir, "gata6_hk.png"),
  gata6_hk_plot,
  height = 12,
  width = 12,
  dpi = 500
)
