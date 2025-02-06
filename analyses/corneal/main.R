## --- CODE BLOCK ---

R_SRC <- Sys.getenv("R_SRC")
M <- new.env()
U <- new.env()
S <- new.env()
P <- new.env()
source(paste0(R_SRC, "/", "utils.R"), local = U)
source(paste0(R_SRC, "/", "sc_rnaseq.R"), local = S)
source(paste0(R_SRC, "/", "plotting.R"), local = P)
source(here::here("analyses", "corneal", ".env.R"), local = M)

ad <- reticulate::import("anndata")
sc_pp <- reticulate::import("scanpy.pp")
sc <- reticulate::import("scanpy")

save_fn <- function(plot, name) {
  ggsave(here(M$outdir, name), plot = plot, dpi = 500, width = 10, height = 10)
}


## --- CODE BLOCK ---
hg37_db <- EnsDb(M$hg37)

qc_helper <- function(adata, db, name, is_adata = TRUE) {
  if (is_adata) {
    sce <- AnnData2SCE(adata)
  } else {
    sce <- adata
  }

  sce <- sce |>
    S$add_feature_info(db, keytype = "GENENAME") |>
    U$x2counts()

  sce <- scuttle::addPerCellQC(sce, subsets = list(mito = rowData(sce)$is_mito))

  mads <- S$qc_mads(sce)
  colData(sce) <- cbind(colData(sce), mads$df)
  sce <- S$identify_reasons(sce, mads$qc_applied)

  plot <- S$plot_qc(sce, c("sum", "detected", "subsets_mito_percent"),
    x_axis = "sample",
    discard_col = "discard_reason"
  )
  save_fn(plot, glue("{name}_qc.png"))
  sce
}


get_fibro <- function(f) {
  adata <- sc$read_10x_mtx(path = here(M$fibro_dir, "Healthy"))
  adata$obs$source <- "fibroblast"
  adata$obs$sample <- "fibroblast_unknown"

  hg37_db <- EnsDb(M$hg37)
  sce <- qc_helper(adata, hg37_db, "fibro")
  writeH5AD(sce, f)
}

get_corneal <- function(f) {
  read <- function(p) {
    adata <- sc$read_10x_mtx(path = here(M$corneal_dir), prefix = p)
    adata$obs$sample <- str_remove(p, "_.*")
    adata
  }
  prefixes <- c("GSM7111220_g016_", "GSM7111221_g017_", "GSM7111223_g019_")
  adata <- ad$concat(lapply(prefixes, read), merge = "same", index_unique = "_")
  adata$obs$source <- "corneal"

  sce <- qc_helper(adata, EnsDb(M$hg38), "corneal")
  writeH5AD(sce, f)
}

get_colon <- function(f) {
  healthy <- c("-A1", "-C1", "-B1")
  df <- read_tsv(here(M$colon_dir, "GSE116222_Expression_matrix.txt"), col_names = FALSE) |> column_to_rownames("X1")
  colnames(df) <- df[1, ]
  df <- df[-1, ]
  wanted_cols <- grepl(paste(healthy, collapse = "|"), colnames(df))
  df <- df[, wanted_cols]
  df <- df |> mutate(across(everything(), as.numeric))

  sce <- SingleCellExperiment::SingleCellExperiment(list(X = df))
  colData(sce)$source <- "colon"
  colData(sce)$sample <- str_remove(colnames(sce), ".*-")
  sce <- qc_helper(sce, EnsDb(M$hg37), "colon", FALSE)
  writeH5AD(sce, f)
}


get_combined <- function(outfile) {
  py$combined <- ad$concat(c(fibro, corneal, colon), merge = "same", join = "outer")
  n_mito <- py$combined$var$is_mito |> sum()

  sc_pp$pca(py$combined)
  sc$external$pp$scanorama_integrate(py$combined, "source")

  sce <- AnnData2SCE(py$combined)
  sce <- scuttle::addPerCellQC(sce,
    subsets = list(mito = rowData(sce)$is_mito),
    assay.type = "X"
  ) |> U$x2counts()

  sce <- scDblFinder(sce, samples = sce$sample)
  writeH5AD(sce, outfile)
  sce
}

# Perform QC individually for each file
fibro_file <- here(M$outdir, "fibro.h5ad")
corneal_file <- here(M$outdir, "corneal.h5ad")
colon_file <- here(M$outdir, "colon.h5ad")

corneal <- U$read_existing(corneal_file, get_corneal, readH5AD)
colon <- U$read_existing(colon_file, get_colon, readH5AD)
fibro <- U$read_existing(fibro_file, get_fibro, readH5AD)

## * Get combined file

## combined_file <- here(M$outdir, "combined_w_ann.h5ad")
## sce <- U$read_existing(combined_file, get_combined, \(x) readH5AD(x, X_name = "counts"))

## # Anotate based on LGR5 positive or negative

## ggplot(
##   as_tibble(colData(sce)),
##   aes(y = subsets_mito_percent, x = source, colour = discard)
## ) +
##   geom_boxplot()

## # Will carry out DGE as pseudobulks with edgeR
## pbs <- Seurat::as.Seurat(sce, data = NULL) |>
##   SeuratObject::RenameAssays("originalexp", "RNA") |>
##   Seurat2PB(sample = "sample", cluster = "source")

## # TODO: build the design matrix

## plot <- P$pca_dgelist(pbs, plot_aes = list(color = "cluster"))
