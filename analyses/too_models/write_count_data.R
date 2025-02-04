source(here::here("analyses", "too_models", "main.R"))

# Write out the count data into gene x sample format with different normalization modes
# and naming schemes

prepare_counts <- function(dge, type) {
  # dge$counts is a gene x sample matrix
  if (type == "edgeR_cpm") {
    # Get cpm without normalization
    cur_counts <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  } else if (type == "edgeR_normalized_cpm") {
    # Get cpm with normalization, within the given tumor type i.e. batch
    cur_counts <- edgeR::cpm(edgeR::normLibSizes(dge), log = TRUE, prior.count = 1)
  } else if (str_detect(type, "fpkm")) {
    gene_lengths <- AnnotationDbi::mapIds(M$db,
      keys = dge$genes[, 1],
      column = "SEQLENGTH", keytype = "GENEID"
    )
    if (type == "log2_fpkm") log <- TRUE
    cur_counts <- DGEobj.utils::convertCounts(dge$counts,
      unit = "fpkm", geneLength = gene_lengths, log = log
    )
  } else if (type == "log2_scaled") {
    cur_counts <- scale(log2(dge$counts + 1), center = TRUE, scale = TRUE)
  } else if (type == "scaled") {
    cur_counts <- scale(dge$counts, center = TRUE, scale = TRUE)
  } else if (type == "log2") {
    cur_counts <- log2(dge$counts)
  } else {
    errorCondition("Type not recognized!")
  }
  cur_counts |>
    as.data.frame() |>
    mutate(gene_id = dge$genes[, 1]) |>
    relocate(gene_id, .before = everything()) |>
    as_tibble()
}

replace_ensembl_ids <- function(df, new_id_col) {
  kept <- colnames(df)
  merged <- inner_join(df, M$id_mapping, by = join_by(x$gene_id == y$ensembl)) |>
    select(-gene_id) |>
    rename(c("gene_id" = new_id_col)) |>
    group_by(gene_id) |>
    mutate(across(where(is.numeric), mean)) |>
    ungroup() |>
    distinct(gene_id, .keep_all = TRUE)
  merged |> select(all_of(kept))
}

QUANT_TYPES <- c(
  ## "edgeR_normalized_cpm",
  ## "kallisto_tpm",
  ## "edgeR_cpm",
  "log2_fpkm",
  "log2_scaled",
  "cupai_log2"
)

unique_tumor_types <- unique(chula_meta$tumor_type)

NAMING_SCHEMES <- list(
  entrez = \(x) replace_ensembl_ids(x, "entrez"),
  ensembl = \(x) x,
  hugo = \(x) replace_ensembl_ids(x, "gene_name")
)

filter_fn <- function(dge) {
  to_keep <- edgeR::filterByExpr(dge, min.count = 10)
  dge[to_keep, ]
}


#' Produce a sample x genes csv, with a column for tumor_type
#'
main_writer <- function(dge, meta, outdir, quant_types, prefix) {
  unique_tumor_types <- unique(meta$tumor_type)
  for (n in names(NAMING_SCHEMES)) {
    for (q in quant_types) {
      all <- list()
      tumor_type <- c()
      dir.create(here(outdir, q))
      for (i in seq_along(unique_tumor_types)) {
        # Separates dataset into tumor types before normalization (relevant to advanced
        # normalization e.g. edgeR)
        t <- unique_tumor_types[i]
        cur_type <- meta |> filter(tumor_type == t)
        cur <- dge[, cur_type$cases] |> filter_fn()
        cur_counts <- prepare_counts(cur, q)
        all[[t]] <- cur_counts
        tumor_type <- c(tumor_type, rep(t, ncol(cur_counts)))
      }
      all_samples <- reduce(all, \(x, y) left_join(x, y, by = join_by(gene_id)))
      to_write <- NAMING_SCHEMES[[n]](all_samples)
      outfile <- here(outdir, q, glue("{prefix}-{q}-{n}.csv"))
      transposed <- U$transpose(to_write, "gene_id") |> cbind(tumor_type)
      write.csv(transposed, outfile)
    }
  }
}

## main_writer(chula_tpm, chula_meta, here(M$out, "kallisto_tpm"), "log2", "chula")
## main_writer(chula_tpm_counts, chula_meta, M$out, c("log2_scaled", "scaled"), "chula")
## main_writer(chula_counts, chula_meta, M$out, c("log2_fpkm", "fpkm"), "chula")

main_writer(tcga_counts, tcga_meta, M$out, c("log2_scaled", "scaled"), "tcga")


## for (n in names(NAMING_SCHEMES)) {
##   for (q in QUANT_TYPES) {
##     all <- list()
##     tumor_type <- c()
##     dir.create(here(M$out, q))
##     for (i in seq_along(unique_tumor_types)) {
##       # Separates dataset into tumor types before normalization (relevant to advanced
##       # normalization e.g. edgeR)
##       t <- unique_tumor_types[i]
##       cur_type <- chula_meta |> filter(tumor_type == t)
##       if (q == "kallisto_tpm") {
##         cur <- chula_tpm[, cur_type$cases] |> filter_fn()
##         cur_counts <- prepare_counts(cur, "log2")
##       } else if (q == "cupai_log2") {
##         cur <- chula_tpm_counts[, cur_type$cases] |> filter_fn()
##         cur_counts <- prepare_counts(cur, "log2_scaled")
##       } else {
##         cur <- chula_counts[, cur_type$cases] |> filter_fn()
##         cur_counts <- prepare_counts(cur, q)
##       }
##       all[[t]] <- cur_counts
##       tumor_type <- c(tumor_type, rep(t, ncol(cur_counts)))
##     }
##     all_samples <- reduce(all, \(x, y) left_join(x, y, by = join_by(gene_id)))
##     to_write <- NAMING_SCHEMES[[n]](all_samples)
##     outfile <- here(M$out, q, glue("chula-{q}-{n}.csv"))
##     transposed <- U$transpose(to_write, "gene_id") |> cbind(tumor_type)
##     write.csv(transposed, outfile)
##   }
## }
