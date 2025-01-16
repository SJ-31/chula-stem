here::i_am("analyses/hcc/hccn_abundance.R")
source(here("analyses", "hcc", "main.R"))
# Attempt to determine fold changes by comparing gene expression between samples directly
# 1. Visualize and correct for batch effects (due to use of public data)
# 2. Carry out between-sample normalization on raw counts
# 3. Compare CPM of MET in P17 and the average of the public normal samples


## * Batch effect visualization

to_keep <- filterByExpr(y = counts)
counts <- counts[to_keep, , keep.lib.sizes = FALSE]

to_keep <- filterByExpr(y = counts)
normalized <- edgeR::normLibSizes(counts)

norm_plot <- pca_dgelist(normalized, plot_aes = list(shape = "group", color = "type"))
raw_plot <- pca_dgelist(counts, plot_aes = list(shape = "group", color = "type"))
save_fn(norm_plot, "norm_pca.png")
save_fn(raw_plot, "raw_pca.png")
# <2025-01-16 Thu> Batch effect is clearly visible

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
# <2025-01-16 Thu> Batch effect reduced slightly, but still not as much as the
# normal version of ComBat
# TODO: can you perform DE with that version? Or just look at the corrected values

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
compare_17 <- function(counts_table, file) {
  n <- counts_table[, M$tcgn_normals] |> rowMeans()
  merged <- bind_cols(n, counts_table[, M$wanted_sample], ccounts$genes) |>
    as_tibble() |>
    `colnames<-`(c("TCGA_normal", "P17", "gene_name"))
  write_tsv(merged, file)
}

normalized_counts <- cpm(ccounts$counts) # CPM
compare_17(normalized_counts, here(outdir, "cpm_P17_tcga_normal.tsv"))
compare_17(combat_corrected, here(outdir, "log2cpm_P17_tcga_normal.tsv"))
