R_SRC <- Sys.getenv("R_SRC")
M <- new.env()
U <- new.env()
S <- new.env()
P <- new.env()
source(paste0(R_SRC, "/", "utils.R"), local = U)
source(paste0(R_SRC, "/", "sc_rnaseq.R"), local = S)
source(paste0(R_SRC, "/", "plotting.R"), local = P)
source(here::here("analyses", "corneal", ".env.R"), local = M)

library(edgeR)
library(Seurat)
options(Seurat.object.assay.version = "v5")

if (path.expand("~") == "/home/shannc") {
  reticulate::use_condaenv("stem-base")
} else {
  reticulate::use_condaenv("nf")
}

ad <- reticulate::import("anndata")
sc_pp <- reticulate::import("scanpy.pp")
sc <- reticulate::import("scanpy")
# %%
CORNEAL_META <- read_csv(here("analyses", "corneal", "metadata.csv"))

## * Helper functions

save_fn <- function(plot, name) {
  ggsave(here(M$outdir, name), plot = plot, dpi = 500, width = 10, height = 10)
}

dimplot_wrapper <- function(obj, file, reduction, group.by) {
  plot <- DimPlot(
    obj,
    reduction = reduction,
    group.by = group.by,
    label = TRUE
  ) &
    NoAxes()
  ggsave(file, plot, width = 20, height = 15, dpi = 500)
}

normalize_to_pca <- function(obj) {
  NormalizeData(
    obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  ) |>
    FindVariableFeatures() |>
    ScaleData() |>
    RunPCA()
}

qc_helper <- function(adata, db, name = NULL, is_adata = TRUE) {
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

  ## sce <- scDblFinder(sce, samples = sce$sample)
  if (!is.null(name)) {
    plot <- S$plot_qc(
      sce,
      c("sum", "detected", "subsets_mito_percent"),
      x_axis = "sample",
      discard_col = "discard_reason"
    )
    save_fn(plot, glue("{name}_qc.png"))
  }

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
  # ";" is the separating character for multiplexed samples
  read <- function(row) {
    p <- row["prefix"]
    adata <- sc$read_10x_mtx(path = here(M$corneal_dir), prefix = p)

    adata$obs$geo <- str_remove(p, "_.*")
    adata$obs$donor <- row["donor"]
    adata$obs$confluency_passage <- row["confluency_passage"]
    adata$obs$library <- row["library"]
    adata$obs$library_effect <- as.character(row["umi_cutoff"])
    adata$obs$day <- row["day"]
    adata$obs$media <- row["media"]

    sce <- qc_helper(adata, EnsDb(M$hg38))
    sce <- sce[, colData(sce)$sum >= row["umi_cutoff"]]
    sce <- sce[, colData(sce)$subsets_mito_percent < 10]
    sce
  }

  sce <- purrr::reduce(apply(CORNEAL_META, 1, read), \(x, y) cbind(x, y))
  # Processing pipeline to recreate results in paper
  obj <- as.Seurat(sce, data = NULL)
  obj <- RenameAssays(obj, "originalexp", "RNA")
  saveRDS(obj, f)
}

get_colon <- function(f) {
  healthy <- c("-A1", "-C1", "-B1")
  ## df <- read_tsv(here("test.tsv"), col_names = FALSE) |> column_to_rownames("X1")
  df <- read_tsv(
    here(M$colon_dir, "GSE116222_Expression_matrix.txt"),
    col_names = FALSE
  ) |>
    column_to_rownames("X1")

  colnames(df) <- df[1, ]
  df <- df[-1, ]
  wanted_cols <- grepl(paste(healthy, collapse = "|"), colnames(df))
  df <- df[, wanted_cols]
  df <- df |>
    mutate(across(everything(), as.numeric)) |>
    as.matrix()

  sce <- SingleCellExperiment::SingleCellExperiment(list(X = df))
  colData(sce)$source <- "colon"
  colData(sce)$sample <- str_remove(colnames(sce), ".*-")
  sce <- qc_helper(sce, EnsDb(M$hg37), "colon", FALSE)

  ## assay(sce, "X") |>
  ##   as.matrix() |>
  ##   Matrix::Matrix(sparse = TRUE)

  writeH5AD(sce, f)
}


get_combined <- function(outfile) {
  py$combined <- ad$concat(
    c(fibro, corneal, colon),
    merge = "same",
    join = "outer"
  )
  n_mito <- py$combined$var$is_mito |> sum()

  sc_pp$pca(py$combined)
  sc$external$pp$scanorama_integrate(py$combined, "source")

  sce <- AnnData2SCE(py$combined)
  sce <- scuttle::addPerCellQC(
    sce,
    subsets = list(mito = rowData(sce)$is_mito),
    assay.type = "X"
  ) |>
    U$x2counts()

  writeH5AD(sce, outfile)
  sce
}

# Perform QC individually for each file
fibro_file <- here(M$outdir, "fibro.h5ad")
corneal_file <- here(M$outdir, "corneal.rds")
colon_file <- here(M$outdir, "colon.h5ad")

corneal <- U$read_existing(corneal_file, get_corneal, readRDS)
## fibro <- U$read_existing(fibro_file, get_fibro, readH5AD)
## colon <- U$read_existing(colon_file, get_colon, readH5AD)

## * Further analysis for corneal

figdir <- here(M$outdir, "corneal_figures")

corneal_integrate <- function(obj) {
  # Get unintegrated clusters for comparison

  # Split by combination of donor and library effect
  obj[[]]$batch <- paste0(obj[[]]$donor, "_", obj[[]]$library_effect)
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)
  obj <- normalize_to_pca(obj)

  obj <- FindNeighbors(obj, reduction = "pca") |>
    FindClusters(algorithm = 1, cluster.name = "clusters_pre") |>
    RunUMAP(dims = 1:30, reduction = "pca", reduction.name = "umap_pre")

  obj <- IntegrateLayers(obj, method = HarmonyIntegration) |>
    FindNeighbors(reduction = "harmony") |>
    FindClusters(algorithm = 1, cluster.name = "clusters_post") |> # 1 is Louvain
    RunUMAP(dims = 1:30, reduction = "harmony", reduction.name = "umap_post")

  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  Idents(obj) <- "clusters_post"
  # They didn't pseudobulk for the DE analysis
  obj
}

integrated_file <- here(M$outdir, "corneal_integrated.rds")
integrated <- U$read_existing(
  integrated_file,
  \(f) {
    obj <- corneal_integrate(corneal)
    saveRDS(obj, f)
    obj
  },
  readRDS
)

markers_v1 <- U$read_existing(
  here(M$outdir, "corneal_markers.csv"),
  \(f) {
    df <- FindAllMarkers(integrated, group.by = "clusters_post")
    write_csv(df, f)
    df
  },
  read_csv
)


## ** Dirty demultiplexing

try_demux <- function(obj, libraries) {
  subset <- obj[, obj[[]]$library %in% libraries]
  subset <- normalize_to_pca(subset)
  obj <- FindNeighbors(obj, reduction = "pca") |>
    FindClusters(algorithm = 1, cluster.name = "batched_clusters") |>
    RunUMAP(dims = 1:30, reduction = "pca", reduction.name = "umap")
  name <- paste0(libraries, collapse = "_")
  dimplot_wrapper(
    obj,
    file = here(outdir, "corneal_demux", glue("{name}_pca.png")),
    reduction = "pca",
    group.by = c("batched_clusters", "donor")
  )
  dimplot_wrapper(
    obj,
    file = here(outdir, "corneal_demux", glue("{name}_umap.png")),
    reduction = "umap",
    group.by = c("batched_clusters", "donor", "library")
  )
  meta <- rownames_to_column(obj[[]], var = "cell") |>
    as_tibble() |>
    dplyr::select(cell, batched_clusters)
  write_csv(meta, file = here(outdir, "corneal_demux", glue("{name}.csv")))
}

try_demux(integrated, c("G1", "G8", "G9"))
try_demux(integrated, c("G2", "G8", "G9"))
try_demux(integrated, c("G3", "G5", "G6", "G7"))
try_demux(integrated, c("G4", "G8", "G9"))


## ** DE analysis

agg <- AggregateExpression(
  integrated,
  group.by = c("clusters_post", "geo"),
  return.seurat = TRUE
)
agg[[]] <- left_join(
  agg[[]],
  select(mutate(CORNEAL_META, geo = str_remove(prefix, "_.*")), -prefix),
  by = join_by(geo)
)

dge <- DGEList(
  counts = as.matrix(LayerData(agg, "counts")),
  samples = agg[[]]
)


## ** Plots

dimplot_wrapper(
  integrated,
  file = here(figdir, "umap_pre.png"),
  reduction = "umap_pre",
  group.by = c("clusters_pre", "donor", "confluency_passage")
)

dimplot_wrapper(
  integrated,
  file = here(figdir, "umap_post.png"),
  reduction = "umap_post",
  group.by = c("clusters_pre", "clusters_post", "donor")
)
