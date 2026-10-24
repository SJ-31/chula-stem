suppressMessages({
  options(
    reticulate.conda_binary = "/home/shannc/Bio_SDD/miniforge3/condabin/conda"
  )
  library(tidyverse)
  library(reticulate)
  library(paletteer)
  library(here)
  source(here("src", "R", "plotting.R"))
  use_condaenv("stem-base")
})

workdir <- here("analyses", "pdac_subtyping")
kras_status <- read_csv(here(workdir, "kras_status.csv")) |>
  mutate(status = factor(status))

ad <- import("anndata")
adata <- ad$read_h5ad(here(workdir, "cohort.h5ad"))
