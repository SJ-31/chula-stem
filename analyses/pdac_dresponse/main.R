library(tidyverse)
library(here)
library(readxl)

# [2025-09-09 Tue] Script for formatting data

dir <- here("analyses", "pdac_dresponse")
response_file <- here(dir, "1-AUC_PDAC.xlsx")


tb_transpose <- function(tb, colname) {
  transposed <- t(tb) |> as.data.frame()
  colnames(transposed) <- transposed[1, ]
  transposed <- transposed[-1, ]
  as_tibble(rownames_to_column(transposed, var = colname))
}

# Sensitive defined as 25% lowest AUC, Resistant 25% highest
sensitivity_cutoffs <- function(x) {
  case_when(
    x < quantile(x, 0.25, na.rm = TRUE) ~ "sensitive",
    x >= quantile(x, 0.75, na.rm = TRUE) ~ "resistant",
    .default = "intermediate"
  )
}

responses <- local({
  r1 <- read_excel(response_file, sheet = 1) |>
    rename_with(\(x) {
      str_remove(str_to_lower(str_replace_all(x, "#", "_")), "case *")
    }) |>
    rename(drug = "...1") |>
    tb_transpose("sample") |>
    mutate(across(-sample, as.numeric))
  r2 <- read_tsv(here(dir, "2025-9-8_paclitaxel.tsv")) |>
    mutate(sample = str_replace_all(sample, "\\.", "_")) |>
    filter(!sample %in% r1$sample)
  r3 <- read_tsv(here(dir, "2025-9-8_gemcitabine.tsv")) |>
    mutate(sample = str_replace_all(sample, "\\.", "_")) |>
    filter(!sample %in% r1$sample)
  reduce(list(r1, r2, r3), bind_rows) |>
    group_by(sample) |>
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) |>
    separate_wider_delim(
      "sample",
      delim = "_",
      names = c("sample", "biopsy"),
      too_few = "align_start"
    )
})
responses$sample <- paste0(
  responses$sample,
  "_",
  replace_na(as.character(responses$biopsy), "")
) |>
  str_remove("_$")
# Manual naming

responses$sample <- case_when(
  str_detect(responses$sample, "PHcase") ~
    str_replace(responses$sample, "PHcase", "PHcase_"),
  str_detect(responses$sample, "^[0-9]+") ~ paste0("P", responses$sample),
  .default = responses$sample
)


li_data <- read_tsv(here(dir, "li_et_al_organoids.tsv"))
li_data <- li_data |>
  filter(Drug.Name %in% colnames(responses)) |>
  tb_transpose("sample") |>
  mutate(across(-sample, as.numeric))
write_tsv(li_data, here(dir, "li_et_al_filtered.tsv"))
write_tsv(responses, here(dir, "all_responses.tsv"))
