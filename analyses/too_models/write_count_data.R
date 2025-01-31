here::i_am("analyses/too_models/write_count_data.R")
source(here::here("analyses", "too_models", "main.R"))

# Write out the count data into gene x sample format with different normalization modes
# and naming schemes


prepare_counts <- function(dge, type) {
  if (type == "edgeR_cpm") {
    # Get cpm without normalization
    cur_counts <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  } else if (type == "edgeR_normalized_cpm") {
    # Get cpm with normalization, within the given tumor type i.e. batch
    cur_counts <- edgeR::cpm(edgeR::normLibSizes(dge),
      log = TRUE, prior.count = 1
    )
  } else if (type == "log2_fpkm") {
    gene_lengths <- AnnotationDbi::mapIds(M$db,
      keys = dge$genes[, 1],
      column = "SEQLENGTH", keytype = "GENEID"
    )
    cur_counts <- DGEobj.utils::convertCounts(dge$counts,
      unit = "fpkm", geneLength = gene_lengths, log = TRUE
    )
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

# <2025-01-30 Thu> Maybe add in normalized counts
quant_types <- c("edgeR_normalized_cpm", "kallisto_tpm", "edgeR_cpm", "log2_fpkm")
unique_tumor_types <- unique(chula_meta$tumor_type)
naming_schemes <- list(
  entrez = \(x) replace_ensembl_ids(x, "entrez"),
  ensembl = \(x) x,
  hugo = \(x) replace_ensembl_ids(x, "gene_name")
)

# Produce a sample x genes csv, with a column for tumor_type
for (n in names(naming_schemes)) {
  for (q in quant_types) {
    all <- list()
    tumor_type <- c()
    dir.create(here(M$out, q))
    for (i in seq_along(unique_tumor_types)) {
      t <- unique_tumor_types[i]
      cur_type <- chula_meta |> filter(tumor_type == t)
      if (q == "kallisto_tpm") {
        cur_counts <- chula_tpm |>
          select(gene_id, all_of(cur_type$cases)) |>
          mutate(across(where(is.numeric), \(x) log2(x + 1)))
      } else {
        cur <- chula_counts[, cur_type$cases]
        cur_counts <- prepare_counts(cur, q)
      }
      all[[t]] <- cur_counts
      tumor_type <- c(tumor_type, rep(t, ncol(cur_counts)))
    }
    all_samples <- reduce(all, \(x, y) left_join(x, y, by = join_by(gene_id)))
    to_write <- naming_schemes[[n]](all_samples)
    outfile <- here(M$out, q, glue("chula-{q}-{n}.csv"))
    transposed <- U$transpose(to_write, "gene_id") |> cbind(tumor_type)
    write.csv(transposed, outfile)
  }
}


## Doing it this way writes the types out individually
## for (i in seq_along(unique_tumor_types)) {
##   t <- unique_tumor_types[i]
##   cur_type <- chula_meta |> filter(tumor_type == t)
##   for (q in quant_types) {
##     dir.create(here(M$out, q))
##     if (q == "kallisto_tpm") {
##       cur_counts <- chula_tpm |>
##         select(gene_id, all_of(cur_type$cases)) |>
##         mutate(across(where(is.numeric), log2))
##     } else {
##       cur <- chula_counts[, cur_type$cases]
##       cur_counts <- prepare_counts(cur, q)
##     }
##     for (n in names(naming_schemes)) {
##       outfile <- here(M$out, q, glue("chula-{q}-{t}-{n}.csv"))
##       to_write <- naming_schemes[[n]](cur_counts)
##       transposed <- U$transpose(to_write, "gene_id") |> cbind(data.frame(tumor_type = t))
##       write.csv(transposed, outfile)
##     }
##   }
## }
