# [2026-02-14 Sat]
# Script to generate adata objects with SCTransform-corrected
# counts, hopefully replicating that what was done to generate the normal samples

# It's probably SCTransform from the filenames, and the fact that SCTransform does
# feature selection, reducing the number of available features. You checked
# and confirmed that the normal files do indeed have fewer features than the raw read in with scanpy
options(future.globals.maxSize = 891289600)

library(tidyverse)
library(reticulate)
library(glue)
library(Seurat)
library(here)
use_condaenv("stem-base")
ad <- import("anndata")
sc <- import("scanpy")

gene_ref <- read_csv(here("analyses", "data", "ensembl_gene_data.csv"))
manifest <- read_csv(here("analyses", "scrnaseq_hcc", "manifest.csv"))
workdir <- here("analyses", "scrnaseq_hcc")
yte <- import("yte")

env <- yte$process_yaml(paste0(
  read_lines(here(workdir, "cellranger_config.yaml")),
  collapse = "\n"
))

# [2026-02-16 Mon] tried version 2 (v2), it didn't look similar
sct_version <- "v1"

sc_transform_to_regress <- c("S.Score", "G2M.Score", "percent.mt")

convert_to_adata <- function(seurat_obj, cols_keep) {
  var <- seurat_obj[["RNA"]][[]][Features(seurat_obj[["SCT"]]), cols_keep]
  ad$AnnData(X = t(seurat_obj[["SCT"]]$counts), var = var)
}

load_and_process <- function(h5_file_path) {
  adata <- sc$read_10x_h5(h5_file_path)
  var_cols_keep <- colnames(adata$var)

  adata$var_names_make_unique()
  adata$var <- left_join(
    adata$var,
    gene_ref,
    by = join_by(x$gene_ids == y$ensembl_gene_id)
  ) |>
    distinct(gene_ids, .keep_all = TRUE)

  adata <- adata[
    ,
    !is.na(adata$var$hgnc_symbol) & !duplicated(adata$var$hgnc_symbol)
  ]
  rownames(adata$var) <- adata$var$hgnc_symbol

  counts <- t(adata$X) |> as.matrix()

  rownames(counts) <- rownames(adata$var)
  obj <- CreateSeuratObject(counts = Matrix::Matrix(counts))
  obj[["RNA"]][[]] <- adata$var

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst")
  obj <- ScaleData(obj, features = rownames(obj))

  obj <- CellCycleScoring(
    obj,
    s.features = cc.genes$s.genes,
    g2m.features = cc.genes$g2m.genes,
    set.ident = TRUE
  )

  obj <- SCTransform(
    obj,
    vars.to.regress = sc_transform_to_regress,
    vst.flavor = sct_version
  )
  convert_to_adata(obj, cols_keep = var_cols_keep)
}

## test_file <- "/home/shannc/Bio_SDD/stem_synology/chula_mount/shannc/output/HCC/SCRNASEQ/HCCN17/processed/cellranger_tumor-control/filtered_feature_bc_matrix.h5"

## obj_test <- load_and_process(test_file)

cr_target <- "filtered_feature_bc_matrix.h5"

. <- apply(manifest, 1, \(row) {
  sn <- row["sample_name"]
  suffix <- glue("{row['type']}-{row['treatment']}")
  outdir <- glue("{env$data_root}/{sn}/processed")
  h5_file <- glue("{outdir}/cellranger_{suffix}/{cr_target}")
  save_file <- glue("{outdir}/{sn}_{suffix}_SCT_{sct_version}.h5ad")
  adata <- load_and_process(h5_file)
  adata$write_h5ad(save_file)
})
