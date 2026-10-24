library(tidyverse)
library(reticulate)
library(pointblank)
library(checkmate)
library(glue)
library(zellkonverter)
library(paletteer)
library(here)
library(SingleCellExperiment)
library(ggplot2)
use_condaenv("stem-base")
source(here("src", "R", "utils.R"))

source(here("src", "R", "plotting.R"))

workdir <- here("analyses", "scrnaseq_hcc")
fn <- new.env()
source_python(here(workdir, "functions.py"), envir = fn)
ad <- import("anndata")
sc <- import("scanpy")
dc <- import("decoupler")
sc_utils <- import("chula_stem.sc_rnaseq")
data <- here("analyses", "data")
# %%

sys <- import("sys")
sys$argv <- c(sys$argv, "test=True")
yte <- import("yte")
env <- yte$process_yaml(paste0(
  read_lines(here(workdir, "cellranger_config.yaml")),
  collapse = "\n"
))
# %%
adata <- ad$read_h5ad(here(workdir, "test", "filtered.h5ad"))

clusters <- ad$read_h5ad(here(
  workdir,
  "test",
  "seurat",
  "harmony_clustering.h5ad"
))

gcol <- "leiden_res1"
tmp <- adata$copy()
tmp$obs <- merge(tmp$obs, clusters$obs, by = "row.names", all.x = TRUE)
tmp <- tmp[!is.na(tmp$obs[[gcol]]), ]

bulked <- dc$pp$pseudobulk(tmp, "sample", groups_col = NULL)
gene_mask <- sc$pp$filter_genes(bulked, min_counts = 50, inplace = FALSE)
cell_mask <- sc$pp$filter_cells(bulked, min_counts = 50, inplace = FALSE)
bulked <- bulked[cell_mask[[1]], gene_mask[[1]]]


library(polars)
library(tidypolars)
library(tidyverse)

scanned <- pl$scan_csv(
  "/home/shannc/Bio_SDD/chula-stem/output/HCC/SCRNASEQ/cohort/annotations/gene_set_activity.csv"
)

val <- scanned$unique(c("group", "clustering"))$group_by(
  "clustering"
)$agg(!!!list(pl$col("group")))$collect()
