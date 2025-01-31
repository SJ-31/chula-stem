library(tximport)

get_meta <- function(workflow_output, type = "STAR") {
  if (type == "STAR") {
    filename <- "_tumor-STAR_Counts.tsv"
  } else if (type == "kallisto") {
    filename <- "_tumor-kallisto.tsv"
  }
  tibble(
    cases = list.files(workflow_output, pattern = "P.*"),
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

workdirs <- c(here(M$remote, "output", "HCC", "RNASEQ"))
## workdirs <- c(
##   here(M$remote, "output", "HCC", "RNASEQ"),
##   here(M$remote, "output", "HCC", "RNASEQ"),
##   here(M$remote, "output", "HCC", "RNASEQ"),
## )
tumor_types <- c("HCC")
## types <- c("HCC", "CCA", "CRC")
metadata <- list()
counts <- list()
tpms <- list()
for (i in seq_along(workdirs)) {
  workdir <- workdirs[i]
  meta <- get_meta(workdir) |> mutate(tumor_type = tumor_types[i])
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
  c <- U$get_rnaseq_counts(meta)

  tpms[[i]] <- tpm
  counts[[i]] <- c
  metadata[[i]] <- meta
}
metadata <- bind_rows(metadata)

counts <- reduce(counts, \(x) left_join(x, y, by = join_by(gene_id)))
tpms <- reduce(tpms, \(x) left_join(x, y, by = join_by(gene_id)))

write_tsv(metadata, M$chula_meta_file)
write_rds(counts, M$chula_raw_counts_file)
write_rds(tpms, M$chula_tpm_file)
