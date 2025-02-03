R_SRC <- Sys.getenv("R_SRC")
U <- new.env()
S <- new.env()
P <- new.env()
source(paste0(R_SRC, "/", "utils.R"), local = U)
source(paste0(R_SRC, "/", "sc_rnaseq.R"), local = S)
source(paste0(R_SRC, "/", "plotting.R"), local = P)
library(tidyverse)
library(here)
library(reticulate)
library(ensembldb)
library(scDblFinder)
library(glue)
library(zellkonverter)
library(sva)

ad <- reticulate::import("anndata")

save_fn <- function(plot, name) {
  ggsave(here(outdir, name), plot = plot, dpi = 500, width = 8, height = 8)
}

M <- list()
M$remote <- here("analyses", "data_all")
M$local <- here("analyses", "data")
M$workflow_output <- here(M$remote, "output", "HCC", "RNASEQ")
M$outdir <- here("analyses", "output", "corneal")

hg38 <- here(M$local, "Homo_sapiens.GRCh38.113.sqlite")
hg37 <- here(M$local, "Homo_sapiens.GRCh37.87.sqlite")


get_combined <- function(outfile) {
  fibro <- readH5AD(fibro_file) |>
    S$add_feature_info(hg37, keytype = "GENENAME") |>
    SCE2AnnData()

  corneal <- readH5AD(corneal_file) |>
    S$add_feature_info(hg37, keytype = "GENENAME") |>
    SCE2AnnData()
  colon <- readH5AD(colon_file) |>
    S$add_feature_info(hg37, keytype = "GENENAME") |>
    SCE2AnnData()
  combined <- ad$concat(c(fibro, corneal, colon), merge = "same")
  print(combined)
  sce <- AnnData2SCE(combined)
  sce <- scuttle::addPerCellQC(sce,
    subsets = list(mito = rowData(sce)$is_mito),
    assay.type = "X"
  ) |> U$x2counts()

  sce <- scDblFinder(sce, samples = sce$sample)

  writeH5AD(sce, outfile)
}

# WARNING: because the cell types are confounded with the batch (i.e. different studies)
# we can't be sure that differences are due to cell type and not batch effect
# TODO: read the experimental protocol to see if they differ strongly with each other

fibro_file <- here(M$outdir, "fibro.h5ad") # GRCh37
corneal_file <- here(M$outdir, "corneal.h5ad") # GRCh38
colon_file <- here(M$outdir, "colon.h5ad") # GRCh37

combined_file <- here(M$outdir, "combined_w_ann.h5ad")
sce <- U$read_existing(combined_file, get_combined, \(x) readH5AD(x, X_name = "counts"))

# TODO: check if the qc thresholds make sense for these cell types
# TODO: check if P' Ta wants you to filter by cell types too
mads <- S$qc_mads(sce, batch_col = "source")
colData(sce) <- cbind(colData(sce), mads$df)

# TODO: make some plots to check the qc filtering is appropriate

# Will carry out DGE as pseudobulks with edgeR
pbs <- Seurat::as.Seurat(sce, data = NULL) |>
  SeuratObject::RenameAssays("originalexp", "RNA") |>
  Seurat2PB(sample = "sample", cluster = "source")

# TODO: build the design matrix

plot <- P$pca_dgelist(pbs, plot_aes = list(color = "cluster"))
