library(tximport)
EXCLUDE <- c("log")

get_meta <- function(workflow_output, type = "STAR") {
  if (type == "STAR") {
    filename <- "_tumor-STAR_Counts.tsv"
  } else if (type == "kallisto") {
    filename <- "_tumor-kallisto.tsv"
  }
  cases <- list.dirs(workflow_output, recursive = FALSE, full.names = FALSE) |>
    discard(\(x) x %in% EXCLUDE) |>
    discard(\(x) str_detect(x, "MultiQC_data"))
  tibble(
    cases = cases,
    files = paste0(
      workflow_output, "/", cases, "/tumor/2-",
      cases, filename
    ),
    source = "Primary Tumor", strandedness = "reverse",
    sequencer = "NovaSeq6000"
  )
}

tx2gene <- AnnotationDbi::select(M$db,
  columns = c("TXID", "GENEID"), keytype = "GENEID",
  keys = keys(M$db)
)

workdirs <- c(
  here(M$remote, "output", "HCC", "RNASEQ"),
  here(M$remote, "output", "CCA", "RNASEQ"),
  here(M$remote, "output", "CRC", "RNASEQ")
)
tumor_types <- c(
  "LIHC",
  "CHOL",
  "COADREAD"
)
metadata <- list()
counts <- list()
tpm_counts <- list()
tpms <- list()

for (i in seq_along(workdirs)) {
  workdir <- workdirs[i]
  meta <- get_meta(workdir, type = "STAR") |> mutate(tumor_type = tumor_types[i])
  k_meta <- get_meta(workdir, type = "kallisto") |> mutate(tumor_type = tumor_types[i])

  tpm <- U$get_rnaseq_counts(k_meta,
    count_col = 2,
    read_fn = \(x) {
      t <- tximport(x,
        tx2gene = tx2gene, importer = read_tsv,
        type = "kallisto", ignoreTxVersion = TRUE
      )
      rownames_to_column(as.data.frame(t$abundance), var = "gene_id")
    }
  )

  counts_from_tpm <- U$get_rnaseq_counts(k_meta, count_col = 2, read_fn = \(x) {
    t <- tximport(x,
      tx2gene = tx2gene, importer = read_tsv,
      type = "kallisto", ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM"
    )
    rownames_to_column(as.data.frame(t$counts), var = "gene_id")
  }) # <2025-02-03 Mon> This is what the CUP-Ai-Dx authors did for their experimental
  # data
  c <- U$get_rnaseq_counts(meta)

  tpms[[i]] <- tpm
  counts[[i]] <- c
  tpm_counts[[i]] <- counts_from_tpm
  metadata[[i]] <- meta
}

metadata <- bind_rows(metadata)
counts <- purrr::reduce(counts, \(x, y) left_join(x, y, by = join_by(gene_id)))
tpms <- purrr::reduce(tpms, \(x, y) left_join(x, y, by = join_by(gene_id)))
tpm_counts <- purrr::reduce(tpm_counts, \(x, y) left_join(x, y, by = join_by(gene_id)))

write_tsv(metadata, M$chula_meta_file)
write_rds(counts, M$chula_raw_counts_file)
write_rds(tpms, M$chula_tpm_file)
write_rds(tpm_counts, M$chula_count_tpm_file)
