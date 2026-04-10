library(tidyverse)
library(glue)
library(here)
library(reticulate)

use_condaenv("stem-base")

data_dir <- "/home/shannc/Bio_SDD/stem_synology/chula_mount/shannc/output/PHCase"
workdir <- here("analyses", "pdac_subtyping")

sample_dirs <- discard(list.files(data_dir), \(x) !str_starts(x, "PHcase_"))

get_purecn_dir <- function(sample_name, root) {
  glue("{root}/{sample_name}/annotations/5-{sample_name}-Call")
}

ad <- import("anndata")

adata <- ad$read_h5ad(here(workdir, "pdac_cohort.h5ad"))

adata$uns$purecn_summary |> write_csv(here(workdir, "purecn_summary.csv"))


get_combined <- function(samples, data_dir, suffix = NULL, log = FALSE) {
  lst <- lapply(samples, \(x) {
    pdir <- get_purecn_dir(x, data_dir)
    if (!is.null(suffix)) {
      read_csv(glue("{pdir}/{x}{suffix}")) |>
        mutate(sample = x) |>
        mutate(across(starts_with("chr"), as.character))
    } else if (log) {
      read_lines(glue("{pdir}/{x}.log"))
    }
  })
  if (!is.null(suffix)) {
    bind_rows(lst)
  } else {
    lst |> `names<-`(samples)
  }
}


all_logs <- get_combined(sample_dirs, data_dir, NULL, TRUE)

warnings <- lapply(all_logs, \(x) {
  x[str_starts(x, "WARN")]
})

unreliable <- lapply(warnings, \(x) {
  x[str_detect(
    x,
    "No copy number alterations found, purity estimate unreliable"
  )]
}) |>
  discard(\(x) length(x) == 0) |>
  names()

## gtb <- get_combined(sample_dirs, data_dir, "_genes.csv") |>
##   filter(!is.na(gene.symbol))
## vtb <- get_combined(sample_dirs, data_dir, "_variants.csv") |>
##   filter(!is.na(gene.symbol))

# Check consistency of ML.C calls
## v_sum <- vtb |>
##   group_by(sample, gene.symbol) |>
##   summarise(n = n(), ML.C = list(unique(ML.C))) |>
##   mutate(n_unique_C = map_dbl(ML.C, length))

## ggplot(tb, aes(x = ML.C, y = C)) +
##   geom_point()
# Compare the maximum likelihood copy number in predictSomatic against that of
# callAlterations

# No predictSomatic output
## tb |>
##   group_by(gene.symbol, sample) |>
##   summarise(C = list(C), ML.C = list(ML.C))
