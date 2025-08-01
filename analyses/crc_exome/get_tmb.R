library(tidyverse)

py_utils <- new.env()
reticulate::source_python(glue("{snakemake@config$src$py}/utils.py"), py_utils)
tmp <- snakemake@output$temp
outfile <- snakemake@output$final

py_utils$parse_multiqc_vep(snakemake@input[[1]], tmp)

read_tsv(tmp) |>
  filter(key %in% ACCEPTED_CONSEQUENCE) |>
  mutate(
    sample = str_extract(sample, "[87]-(.*)-VEP_small"),
    key = map_chr(key, \(x) str_replace_all(str_to_title(x), "_", " "))
  ) |>
  filter(category == "Consequences (all)") |>
  write_tsv(outfile)
