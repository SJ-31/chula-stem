R_SRC <- Sys.getenv("R_SRC")
M <- new.env()
U <- new.env()
S <- new.env()
P <- new.env()
source(paste0(R_SRC, "/", "utils.R"), local = U)
source(paste0(R_SRC, "/", "sc_rnaseq.R"), local = S)
source(paste0(R_SRC, "/", "plotting.R"), local = P)
source(here::here("analyses", "corneal", ".env.R"), local = M)
figdir <- here(M$outdir, "corneal_figures")

library(edgeR)
library(tidyverse)
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


corneal_file <- here(M$outdir, "corneal.rds")

corneal <- U$read_existing(corneal_file, get_corneal, readRDS)

## * Analyses

## ** Integrate
corneal_integrate <- function(obj) {
  # Get unintegrated clusters for comparison

  # Split by combination of donor and library effect
  obj[[]]$batch <- paste0(obj[[]]$donor, "_", obj[[]]$library_effect)
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)
  obj <- normalize_to_pca(obj)

  obj <- FindNeighbors(obj, reduction = "pca") |>
    FindClusters(algorithm = 1, cluster.name = "clusters_pre") |>
    RunUMAP(dims = 1:30, reduction = "pca", reduction.name = "umap_pre")

  obj <- IntegrateLayers(
    obj,
    method = HarmonyIntegration
  ) |>
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
integrated[[]]$seurat_clusters <- NULL
integrated[[]]$cluster_confluence <- paste0(
  integrated[[]]$clusters_post,
  "_passage_",
  integrated[[]]$confluency_passage
)
integrated <- integrated[, !integrated[[]]$library %in% paste0("G", 1:4)]
# NOTE: the above is temporary, discarding the multiplexed samples until you can annotate them
rm(corneal)

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
  subset <- FindNeighbors(subset, reduction = "pca") |>
    FindClusters(algorithm = 1, cluster.name = "batched_clusters") |>
    RunUMAP(dims = 1:30, reduction = "pca", reduction.name = "umap")
  name <- paste0(libraries, collapse = "_")
  dimplot_wrapper(
    subset,
    file = here(M$outdir, "corneal_demux", glue("{name}_pca.png")),
    reduction = "pca",
    group.by = c("batched_clusters", "donor")
  )
  dimplot_wrapper(
    subset,
    file = here(M$outdir, "corneal_demux", glue("{name}_umap.png")),
    reduction = "umap",
    group.by = c("batched_clusters", "donor", "library")
  )
  meta <- rownames_to_column(subset[[]], var = "cell") |>
    as_tibble() |>
    dplyr::select(cell, batched_clusters)
  write_csv(meta, file = here(M$outdir, "corneal_demux", glue("{name}.csv")))
}

do_demux <- FALSE
if (path.expand("~") != "/home/shannc" && do_demux) {
  try_demux(corneal, c("G1", "G8", "G9"))
  try_demux(corneal, c("G2", "G8", "G9"))
  try_demux(corneal, c("G3", "G5", "G6", "G7"))
  try_demux(corneal, c("G4", "G8", "G9"))
}

## ** DE analysis

agg <- AggregateExpression(
  integrated,
  group.by = c("clusters_post", "geo"),
  return.seurat = TRUE
)
agg[[]] <- left_join(
  agg[[]],
  dplyr::select(mutate(CORNEAL_META, geo = str_remove(prefix, "_.*")), -prefix),
  by = join_by(geo)
)

dge <- DGEList(
  counts = as.matrix(LayerData(agg, "counts")),
  samples = agg[[]]
)


## ** Plots

U$read_existing(here(figdir, "umap_pre.png"), \(f) {
  dimplot_wrapper(
    integrated,
    file = f,
    reduction = "umap_pre",
    group.by = c("clusters_pre", "donor", "confluency_passage")
  )
})

U$read_existing(here(figdir, "umap_post.png"), \(f) {
  dimplot_wrapper(
    integrated,
    file = f,
    reduction = "umap_post",
    group.by = c("clusters_pre", "clusters_post", "donor")
  )
})


## *** Original markers

markers <- list(
  proliferative = c("CENPF", "PTTG1", "MKI67", "PCNA", "TOP2A"),
  senescent = c("MT2A", "CDKN2A", "TAGLN"),
  fibrotic = c("ACTA2", "CD44", "COL6A1", "COL6A3"),
  cec_phenotype = c("COL4A3", "CDH2", "ALCAM", "SLC4A11"),
  ecm_activity = c("COL4A1", "COL4A2", "COL5A1", "FBLN5")
)


U$read_existing(here(figdir, "original_markers_cluster_confluency.png"), \(f) {
  plot <- DotPlot(
    integrated,
    markers,
    group.by = "cluster_confluence",
    cols = c("lightgrey", "blue")
  ) +
    ylab("Cluster, Confluency Passage")
  ggsave(f, plot, width = 20, height = 15, dpi = 500)
})

U$read_existing(here(figdir, "original_markers_cluster.png"), \(f) {
  plot <- DotPlot(
    integrated,
    markers,
    cols = c("lightgrey", "blue")
  ) +
    ylab("Cluster")
  ggsave(f, plot, width = 20, height = 15, dpi = 500)
})

# [2025-07-02 Wed]
# WARNING: the dotplot is misleading: the representation of cluster 12
# for genes PTTG1, CENPF being highly expressed is false
# you confirmed this with the violin plot, average expression and FindAllMarkers results

gene2class <- tibble(gene = markers, class = names(markers)) |>
  unnest(cols = c(gene))
average_expr <- AggregateExpression(
  integrated,
  features = unlist(markers)
)$RNA |>
  as.data.frame() |>
  rownames_to_column(var = "gene") |>
  as_tibble() |>
  left_join(gene2class, by = join_by(gene))

U$read_existing(here(figdir, "original_markers_violin.png"), \(f) {
  plot <- VlnPlot(
    integrated,
    c(markers$proliferative, markers$senescent),
    group.by = "clusters_post",
  )
  ggsave(f, plot, width = 15, height = 15, dpi = 500)
})

## * Proliferative markers
p_clust <- 11 # [2025-07-02 Wed] Cluster containing proliferative cells

# top markers from FindAllMarkers
up_fam <- markers_v1 |>
  dplyr::filter(cluster == p_clust) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  pluck("gene")

VlnPlot(integrated, features = up_fam)


p_conserved <- FindConservedMarkers(
  # Markers conserved across all passages
  integrated,
  11,
  grouping.var = "confluency_passage"
)
clusters <- integrated[[]]$clusters_post

dge <- normLibSizes(dge)
mm <- model.matrix(~ clusters_post + confluency_passage, data = colData(dge))
dge <- estimateDisp(dge, design = mm, robust = TRUE)
fit <- glmQLFit(dge, mm, robust = TRUE)
mean_val <- 1 / (length(clusters) - 1)
mean_others <- paste0(
  mean_val,
  "*",
  clusters[clusters != p_clust],
  collapse = "+"
)
contrast_str <- paste0(clusters, "-", "(", mean_others, ")")
test <- glmQLFTest(fit, contrast = )
