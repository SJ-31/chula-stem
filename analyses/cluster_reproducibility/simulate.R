library(here)
library(glue)
library(tidyverse)
library(SPsimSeq)
library(reticulate)
use_condaenv("stem-base")

yte <- import("yte")
env <- yte$process_yaml(paste0(
  read_lines(here("analyses", "cluster_reproducibility", "env.yaml")),
  collapse = "\n"
))
seed <- env$seed

ad <- import("anndata")

workdir <- do.call(here, as.list(env$workdir))
outdir <- do.call(here, as.list(env$outs$simulated))
dir.create(outdir)


data("zhang.data.sub")
zhang_counts <- zhang.data.sub$counts
group <- zhang.data.sub$MYCN.status

sim_spec <- env$datasets$simulation

for (dset in names(sim_spec)) {
  spec <- sim_spec[[dset]]
  grouping <- sample(
    letters[1:n_groups],
    size = spec$tot.samples,
    replace = TRUE
  )
  sim_result <- do.call(
    \(...) {
      SPsimSeq(
        n.sim = 1,
        s.data = zhang_counts,
        group = group,
        result.format = "list",
        group.config = c(0.5, 0.5),
        ...
      )
    },
    spec
  )
  counts <- sim_result[[1]]$counts
  var <- sim_result[[1]]$rowData
  adata <- ad$AnnData(X = t(counts), var = var)
  adata$write_h5ad(here(outdir, glue("{dset}.h5ad")))
}
