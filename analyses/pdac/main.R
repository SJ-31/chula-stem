library(tidyverse)
library(ggplot2)
library(paletteer)
library(glue)
library(here)
library(reticulate)
library(ggpubr)
library(cowplot)

utils <- new.env()
source(here("src", "R", "utils.R"), local = utils)
py_utils <- new.env()
reticulate::source_python(here("src", "chula_stem", "utils.py"), py_utils)

vaf_merged_file <- here("analyses", "output", "pdac_vaf_merged.tsv")
sbs_merged_file <- here("analyses", "output", "pdac_sbs.tsv")
sbs_merged_file_mp <- here("analyses", "output", "pdac_sbs_mp.tsv")
sbs_merged_file_all <- here("analyses", "output", "pdac_sbs_all.tsv")
vaf_merged_transcripts_file <- here("analyses", "output", "pdac_vaf_merged_transcripts.tsv")
data_path <- here("analyses", "data_all", "output", "PDAC")

# TODO: when pdac is done, recreate this file
if (!file.exists(vaf_merged_file)) {
  files <- list.files(data_path, pattern = "8-P[0-9_]+-VEP_small.tsv$", recursive = TRUE, full.names = TRUE)
  tsvs <- lapply(files, \(x) {
    tb <- read_with_filename(x, "sample") |> select(sample, VAF, Feature, SYMBOL, CLIN_SIG, Consequence, SOURCE)
    tb |>
      mutate(census_mutations = map_chr(Consequence, \(x) {
        case_when(
          str_detect(x, "missense") ~ "Mis",
          str_detect(x, "frameshift") ~ "F",
          str_detect(x, "deletion|lost") ~ "D",
          .default = NA
        )
      })) |>
      group_by(SYMBOL) |>
      summarise(
        VAF = mean(VAF),
        across(where(is.character), \(x) utils$flatten_by(x, collapse = TRUE, unique = FALSE))
      ) |>
      mutate(census_mutations = map_chr(census_mutations, \(x) {
        if (is.na(x)) {
          return(x)
        }
        first(utils$modes(str_split_1(x, ";") |> discard(\(x) x == "NA")))
      }))
  })
  merged <- bind_rows(tsvs) |> mutate(sample = str_extract(sample, "P[0-9_]+"))
  write_tsv(merged, vaf_merged_file)
  q()
} else {
  merged <- read_tsv(vaf_merged_file) |>
    mutate(across(c(Consequence, CLIN_SIG), utils$into_char_list))
}

multiqc_file <- here(data_path, "8-cohort-MultiQC_data", "vep.txt")
vep_data_file <- here("analyses", "output", "pdac_vep_data.tsv")
if (!file.exists(vep_data_file)) {
  py_utils$parse_multiqc_vep(multiqc_file, vep_data_file)
}
vep_data <- read_tsv(vep_data_file)
