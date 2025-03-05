source(here::here("analyses", "too_models", "main.R"))

fpkm_uq <- function(dge, lengths) {
  # Calculation from
  # https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#mrna-expression-transformation
  quantiles <- apply(dge$counts, 2, \(x) quantile(x, 0.75))
  qdf <- matrix(quantiles, nrow = nrow(dge), ncol = ncol(dge), byrow = TRUE)
  n_autosomes <- 19029
  denom <- (qdf * n_autosomes * lengths)
  normalized <- (dge$counts * 10^9) / denom
  dge$counts <- normalized
  dge
}

# Write out the count data into gene x sample format with different normalization modes
# and naming schemes

prepare_counts <- function(dge, type) {
  # dge$counts is a gene x sample matrix
  gene_lengths <- AnnotationDbi::mapIds(M$db,
    keys = dge$genes[, 1],
    column = "SEQLENGTH", keytype = "GENEID"
  )
  if (str_detect(type, "log2_")) {
    log <- TRUE
  } else {
    log <- FALSE
  }

  if (type == "edgeR_cpm") {
    # Get cpm without normalization
    cur_counts <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  } else if (type == "edgeR_normalized_cpm") {
    # Get cpm with normalization, within the given tumor type i.e. batch
    cur_counts <- edgeR::cpm(edgeR::normLibSizes(dge), log = TRUE, prior.count = 1)
  } else if (str_detect(type, "fpkm")) {
    if (type == "fpkm_uq") {
      cur_counts <- fpkm_uq(dge, gene_lengths)
    } else {
      cur_counts <- DGEobj.utils::convertCounts(dge$counts,
        unit = "fpkm", geneLength = gene_lengths, log = log, prior.count = 1
      )
    }
  } else if (type == "log2_scaled") {
    cur_counts <- scale(log2(dge$counts + 1), center = TRUE, scale = TRUE)
  } else if (str_detect(type, "tpm")) {
    cur_counts <- DGEobj.utils::convertCounts(dge$counts,
      unit = "TPM", geneLength = gene_lengths,
      log = log, prior.count = 1
    )
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
  tulip = \(x) {
    tb <- x |>
      dplyr::select(gene_id, where(is.numeric)) |>
      left_join(M$id_mapping, by = join_by(x$gene_id == y$ensembl)) |>
      dplyr::select(-entrez) |>
      relocate(gene_name, .after = gene_id)
    tb
  },
  hugo = \(x) replace_ensembl_ids(x, "gene_name")
)

filter_fn <- function(dge) {
  to_keep <- edgeR::filterByExpr(dge, min.count = 10)
  dge[to_keep, ]
}


#' Produce a sample x genes csv, with a column for tumor_type
#'
main_writer <- function(dge, meta, outdir, quant_types, prefix, filter = TRUE,
                        name_as = NULL, transpose = TRUE) {
  unique_tumor_types <- unique(meta$tumor_type)
  if (is.null(name_as)) {
    name_as <- names(NAMING_SCHEMES)
  }
  for (n in name_as) {
    for (q in quant_types) {
      all <- list()
      tumor_type <- c()
      dir.create(here(outdir, q))
      for (i in seq_along(unique_tumor_types)) {
        # Separates dataset into tumor types before normalization (relevant to advanced
        # normalization e.g. edgeR)
        t <- unique_tumor_types[i]
        cur_type <- meta |> dplyr::filter(tumor_type == t)
        if (filter) {
          cur <- dge[, cur_type$cases] |> filter_fn()
        } else {
          cur <- dge[, cur_type$cases]
        }
        cur_counts <- prepare_counts(cur, q)
        all[[t]] <- cur_counts
        tumor_type <- c(tumor_type, rep(t, ncol(cur_counts)))
      }
      all_samples <- purrr::reduce(all, \(x, y) left_join(x, y, by = join_by(gene_id)))
      to_write <- NULL
      try(to_write <- NAMING_SCHEMES[[n]](all_samples), silent = TRUE)
      outfile <- here(outdir, q, glue("{prefix}-{q}-{n}.csv"))
      if (!is.null(to_write)) {
        warning(glue("Attempt to write {outfile} failed!"))
        if (transpose) {
          transposed <- U$transpose(to_write, "gene_id") |> cbind(tumor_type)
          write.csv(transposed, outfile)
        } else {
          write_csv(to_write, outfile)
        }
      }
    }
  }
}

## chula_tpm <- read_rds(M$chula_tpm_file) |>
##   DGEList(samples = chula_meta, group = chula_meta$tumor_type)
## main_writer(chula_tpm, chula_meta, M$out, "log2", "chula_tpm")

## chula_tpm_counts <- read_rds(M$chula_count_tpm_file) |>
##   DGEList(samples = chula_meta, group = chula_meta$tumor_type)
## main_writer(chula_tpm_counts, chula_meta, M$out, c("log2_scaled", "scaled"), "chula")

## main_writer(chula_counts, chula_meta, M$out, c("log2_fpkm", "fpkm"), "chula")

## main_writer(tcga_counts, tcga_meta, M$out, c(
##   "log2_scaled",
##   "scaled",
##   "log2_tpm"
##   ## "tpm",
##   ## "log2_fpkm"
## ), "tcga")

firehose_counts <- readRDS(here(M$out, "firehose_counts.rds"))
firehose_meta <- as_tibble(firehose_counts$samples)
main_writer(firehose_counts, firehose_meta, M$out, c("log2", "tpm", "log2_fpkm", "log2_tpm"), "firehose")

public_counts <- readRDS(here(M$out, "public_counts.rds"))
public_meta <- as_tibble(public_counts$samples)
main_writer(public_counts, public_meta, M$out, c("log2_tpm", "tpm", "log2_fpkm"), "public")

## main_writer(tcga_counts, tcga_meta, M$out, "fpkm_uq", "tcga", name_as = "tulip", transpose = FALSE)

# <2025-02-04 Tue> Try to quantile-normalize the tcga data, because this appears to be
# what Firehose does
# See https://broadinstitute.atlassian.net/wiki/spaces/GDAC/pages/844334346/Documentation
qn_tcga <- qn_counts(tcga_counts)
main_writer(qn_tcga, tcga_meta, M$out, c(
  "log2_scaled", "scaled"
), "qn-tcga", filter = FALSE)
