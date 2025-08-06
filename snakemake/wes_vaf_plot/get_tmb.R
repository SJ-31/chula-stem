library(tidyverse)
library(glue)
if (!is.null(snakemake@config$conda)) {
  reticulate::use_condaenv(snakemake@config$conda)
}
py_utils <- new.env()
reticulate::source_python(glue("{snakemake@config$src$py}/utils.py"), py_utils)
tmp <- snakemake@output$temp
outfile <- snakemake@output[[1]]

lapply(snakemake@input, \(x) {
  py_utils$parse_multiqc_vep(x) |>
    as_tibble() |>
    mutate(sample = str_extract(sample, "[87]-(.*)-VEP_small", group = 1)) |>
    filter(category == "Consequences (all)")
}) |>
  bind_rows() |>
  write_tsv(outfile)
