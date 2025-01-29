library(tidyverse)
library(ggplot2)
library(paletteer)
library(glue)
library(here)
library(reticulate)
library(cowplot)

here::i_am("./analyses/pdac/main.R")
utils <- new.env()
source(here("src", "R", "utils.R"), local = utils)
py_utils <- new.env()
reticulate::source_python(here("src", "chula_stem", "utils.py"), py_utils)

outdir <- here("analyses", "output", "pdac")
save_fn <- function(plot, name) {
  ggsave(here(outdir, name), plot = plot, dpi = 500, width = 15, height = 10)
}


vaf_merged_file <- here("analyses", "output", "pdac_vaf_merged.tsv")
sbs_merged_file <- here("analyses", "output", "pdac", "sbs.tsv")
sbs_merged_file_mp <- here("analyses", "output", "pdac_sbs_mp.tsv")
sbs_merged_file_all <- here("analyses", "output", "pdac_sbs_all.tsv")
vaf_merged_transcripts_file <- here("analyses", "output", "pdac_vaf_merged_transcripts.tsv")
data_path <- here("analyses", "data_all", "output", "PDAC")
ref_path <- here("analyses", "data_all", "reference")
dbsnp_file <- here(ref_path, "variants", "dbSNP_somatic.csv")

# TODO: when pdac is done, recreate this file
if (!file.exists(vaf_merged_file)) {
  files <- list.files(data_path, pattern = "8-P[0-9_]+-VEP_small.tsv$", recursive = TRUE, full.names = TRUE)
  tsvs <- lapply(files, \(x) {
    tb <- utils$read_with_filename(x, "sample") |> select(
      sample, Alt_depth,
      VAF, Feature, SYMBOL, Existing_variation,
      CLIN_SIG, Consequence, SOURCE
    )
    ## tb |> <2025-01-24 Fri> Takes too long and you don't use this anyway
    ##   mutate(census_mutations = map_chr(Consequence, \(x) {
    ##     case_when(
    ##       str_detect(x, "missense") ~ "Mis",
    ##       str_detect(x, "frameshift") ~ "F",
    ##       str_detect(x, "deletion|lost") ~ "D",
    ##       .default = NA
    ##     )
    ##   })) |>
    ##   mutate(across(where(is.character), \(x) utils$flatten_by(x, collapse = TRUE, unique = FALSE))) |>
    ##   mutate(census_mutations = map_chr(census_mutations, \(x) {
    ##     if (is.na(x)) {
    ##       return(x)
    ##     }
    ##     first(utils$modes(str_split_1(x, ";") |> discard(\(x) x == "NA")))
    ##   }))
    tb
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

filter_known <- function(tb, dbsnp) {
  if (is.character(dbsnp)) {
    dbsnp <- read_csv(dbsnp_file)
  }
  cur_cols <- colnames(tb) |> setdiff("Existing_variation")
  tb |>
    separate_longer_delim(Existing_variation, ";") |>
    filter(Existing_variation %in% dbsnp$ID | grepl("COSV", Existing_variation)) |>
    distinct(pick(all_of(cur_cols)))
}

vaf_heatmap <- function(plot) {
  plot +
    geom_tile(width = 0.95, height = 0.95) +
    xlab("Sample") + ylab("Gene") +
    theme_grey() +
    theme(
      plot.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.ticks = element_line(size = 0.5),
      axis.text.x = element_text(angle = 90),
    )
}
