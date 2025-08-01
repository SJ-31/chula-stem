library(tidyverse)
library(glue)
reticulate::use_condaenv(snakemake@config$conda)
py_utils <- new.env()
reticulate::source_python(glue("{snakemake@config$src$py}/utils.py"), py_utils)
tmp <- snakemake@output$temp
outfile <- snakemake@output[[1]]

df <- py_utils$parse_multiqc_vep(snakemake@input[[1]]) |> as_tibble()


accepted_consequence <- snakemake@config$accepted_consequence
if (is.null(accepted_consequence)) {
  accepted_consequence <- df$key
}

df |>
  filter(key %in% accepted_consequence) |>
  mutate(
    sample = str_extract(sample, "[87]-(.*)-VEP_small"),
    key = map_chr(key, \(x) str_replace_all(str_to_title(x), "_", " "))
  ) |>
  filter(category == "Consequences (all)") |>
  write_tsv(outfile)
