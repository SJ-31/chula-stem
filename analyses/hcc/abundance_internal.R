here::i_am("analyses/hcc/abundance_internal.R")
source(here::here("analyses", "hcc", "main.R"))
library(DGEobj.utils)

# Attempt to determine fold changes between tumor and normal internally, with housekeeping
# genes
# 1. Select appropriate housekeeping genes (HG) for HCC
# 2. Carry out within-sample normalization on raw counts
# 3. For each sample, compare relative log2 fold change in MET expression with HGs
# 4. Compare differences between tumor and normal samples
#   Assumption is that the expression of HGs does not change regardless of condition, so
#   any differences in fold change are due to MET expression change
gene_lengths_file <- here(outdir, "gene_lengths.tsv")
if (!file.exists(gene_lengths_file)) {
  gene_lengths <- local({
    data <- read_rds(here(data_dir, "reference", "genomes", "Homo_sapiens.GRCh38.113.rds"))
    lengths <- as_tibble(data) |>
      select(type, gene_name, start, end) |>
      filter(type == "gene") |>
      distinct(gene_name, .keep_all = TRUE) |>
      mutate(length = end - start) |>
      select(gene_name, length) |>
      filter(!is.na(gene_name))
    lengths
  })
  write_tsv(gene_lengths, gene_lengths_file)
} else {
  gene_lengths <- read_tsv(gene_lengths_file)
}

## * Convert raw counts to TPM
to_keep <- filterByExpr(y = counts, min.count = 10)
counts <- counts[to_keep, , keep.lib.sizes = FALSE]

aligned_length <- left_join(counts$genes, gene_lengths, by = join_by(gene_name))
no_lengths <- filter(aligned_length, is.na(length)) |> pluck("gene_name")
unknown_removed <- counts[!counts$genes[, 1] %in% no_lengths, , keep.lib.sizes = FALSE]
unknown_removed$genes <- left_join(unknown_removed$genes, aligned_length, by = join_by(gene_name))

chula_samples <- counts$samples |>
  filter(group == "chula") |>
  rownames()
samples <- c(M$tcgn_normals, M$wanted_sample)
samples <- c(M$tcgn_normals, chula_samples)
tpm <- DGEobj.utils::convertCounts(unknown_removed$counts, "TPM", geneLength = unknown_removed$genes$length) |>
  as.data.frame() |>
  dplyr::select(all_of(samples)) |>
  bind_cols(unknown_removed$genes) |>
  `rownames<-`(NULL) |>
  column_to_rownames(var = "gene_name")

tpm_not_zero <- tpm |>
  filter(if_all(is.numeric, \(x) x != 0)) |>
  rownames()

## * Compute with-sample fold log2_change against housekeeping genes
fold_change <- function(target, housekeeping, data, samples = NULL,
                        const = 1,
                        verbose = FALSE, log2 = FALSE) {
  # Percent change col means that `target`'s expression is x% of the hg
  helper <- function(h) {
    if (is.null(samples)) {
      s <- colnames(data)
    } else {
      s <- samples
    }
    log2_change <- map_dbl(s, \(x) {
      expr_target <- data[target, x]
      if (verbose) {
        print(glue("\n\n---\t{target} tpm in {x}: {expr_target}"))
      }
      diff <- log(expr_target + const) - log(data[h, x] + const)
      if (verbose) print(glue("{h} tpm: {data[h, x]}"))
      diff
    })
    tb <- tibble(sample = s)
    if (!log2) {
      tb[[h]] <- 2^(log2_change) * 100
    } else {
      tb[[h]] <- log2_change
    }
    tb
  }
  reduce(lapply(housekeeping, helper), \(x, y) inner_join(x, y, by = join_by(sample)))
}

# Top two HGs from the HRT Atlas website https://housekeeping.unicamp.br/?homePageHuman
hg1 <- c("GABARAP", "PFDN5")

# Top two HGs based on https://www.spandidos-publications.com/10.3892/ol.2021.13052#f1-ol-0-0-13052
hg2 <- c("HMBS", "RPLP0")

# Most stable gene from `hcc_abundance.R` script,
# whose tpm are all nonzero
# Using TCGA data as of <2025-01-16 Thu>
hg3 <- c("GNAI3")

all_hg <- c(hg1, hg2, hg3)

samples <- colnames(tpm) |> discard(\(x) x %in% c("length"))
met_change <- fold_change("MET", all_hg, tpm, samples = samples, verbose = TRUE)
met_change_file <- here(outdir, "met_fold_change_internal.tsv")

test_change <- wilcox.test(
  x = filter(met_change, sample == M$wanted_sample)[, 3][[1]],
  y = filter(met_change, sample == M$tcgn_normals)[, 3][[1]],
  data = ""
)

longer <- pivot_longer(met_change, cols = where(is.numeric)) %>%
  inner_join(counts$samples, by = join_by(x$sample == y$files)) |>
  filter(sample == M$wanted_sample | tissue == "N")

fold_plot <- ggplot(longer, aes(y = value, x = name, color = tissue)) +
  geom_boxplot() +
  xlab("Housekeeping gene") +
  ylab("MET fold-change (%)")

write_tsv(met_change, met_change_file)
save_fn(fold_plot, "met_fold_change.png")

# TODO: look at the data TMM normalized data to find housekeeping genes manually
