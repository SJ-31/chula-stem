library(zellkonverter)
library(scRNAseq)
source(here::here("analyses", "too_models", "paths.R"))

data <- readH5AD(here(out, "TCGA_COAD-READ_tpm.h5ad"))

ids <- str_remove_all(rowData(data)$gene_id, "\\..*")
rowData(data)$gene_id <- ids
data <- data[!duplicated(rowData(data)$gene_id), ]

entrez <- left_join(
  data.frame(gene_id = rowData(data)$gene_id),
  id_mapping,
  by = join_by(x$gene_id == y$ensembl),
) |> distinct(gene_id, .keep_all = TRUE)

rowData(data)$entrez <- pluck(entrez, "entrez")
data <- data[!is.na(rowData(data)$entrez), ]

counts <- log(t(assays(data)$X) + 1) |> as.data.frame()
colnames(counts) <- rowData(data)$entrez
rownames(counts) <- colData(data)$File_ID
counts$tumor_type <- str_remove(colData(data)$Project_ID, "TCGA-")

write.csv(counts, here(out, "tcga_coad-read2.csv"))
