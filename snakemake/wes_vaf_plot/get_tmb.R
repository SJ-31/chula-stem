library(tidyverse)
library(glue)
outfile <- snakemake@output[[1]]

parse_multiqc_vep <- function(file) {
  cols <- c(
    "Variant classes",
    "Consequences (most severe)",
    "Consequences (all)",
    "Coding consequences",
    "Variants by chromosome",
    "Position in protein",
    "General statistics"
  )
  read_tsv(file) |>
    mutate(across(all_of(cols), \(col) {
      lapply(col, \(x) {
        rep <- str_replace_all(x, "'", "\"")
        jsonlite::parse_json(rep)
      })
    })) |>
    pivot_longer(-Sample, names_to = "category") |>
    rename(sample = "Sample") |>
    unnest(value) |>
    mutate(key = names(value), value = unlist(value, use.names = FALSE))
}


lapply(snakemake@input, \(x) {
  parse_multiqc_vep(x) |>
    mutate(sample = str_extract(sample, "[87]-(.*)-VEP_small", group = 1)) |>
    filter(category == "Consequences (all)")
}) |>
  bind_rows() |>
  write_tsv(outfile)
