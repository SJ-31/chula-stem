library(tidyverse)
library(glue)
reticulate::use_condaenv(snakemake@config$conda)
py_utils <- new.env()
reticulate::source_python(glue("{snakemake@config$src$py}/utils.py"), py_utils)
tmp <- snakemake@output$temp
outfile <- snakemake@output[[1]]

py_utils$parse_multiqc_vep(snakemake@input[[1]]) |>
  as_tibble() |>
  mutate(sample = str_extract(sample, "[87]-(.*)-VEP_small")) |>
  filter(category == "Consequences (all)") |>
  write_tsv(outfile)
