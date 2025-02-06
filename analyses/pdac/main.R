library(tidyverse)
library(ggplot2)
library(paletteer)
library(glue)
library(here)
library(reticulate)
library(cowplot)

utils <- new.env()
source(here("src", "R", "utils.R"), local = utils)
py_utils <- new.env()
reticulate::source_python(here("src", "chula_stem", "utils.py"), py_utils)

outdir <- here("analyses", "output", "pdac")
save_fn <- function(plot, name) {
  ggsave(here(outdir, name), plot = plot, dpi = 500, width = 15, height = 10)
}

vaf_merged_file <- here(outdir, "pdac_vaf_merged.tsv")
sbs_merged_file <- here(outdir, "sbs.tsv")

vaf_merged_transcripts_file <- here(outdir, "pdac_vaf_merged_transcripts.tsv")
data_path <- here("analyses", "data_all", "output", "PDAC")
ref_path <- here("analyses", "data_all", "reference")
dbsnp_file <- here("analyses", "data", "dbSNP_somatic.csv")

if (!file.exists(vaf_merged_file)) {
  files <- list.files(data_path, pattern = "8-P[0-9_]+-VEP_small.tsv$", recursive = TRUE, full.names = TRUE)
  tsvs <- lapply(files, \(x) {
    tb <- utils$read_with_filename(x, "sample") |> select(
      sample, Alt_depth,
      VAF, Feature, SYMBOL, Existing_variation,
      CLIN_SIG, Consequence, SOURCE
    )
    tb
  })
  merged <- bind_rows(tsvs) |> mutate(sample = str_extract(sample, "P[0-9_]+"))
  write_tsv(merged, vaf_merged_file)
} else {
  merged <- read_tsv(vaf_merged_file) |>
    mutate(across(c(Consequence, CLIN_SIG), utils$into_char_list))
}

multiqc_file <- here(outdir, "vep.txt")
vep_data_file <- here(outdir, "pdac_vep_data.tsv")
if (!file.exists(vep_data_file)) {
  py_utils$parse_multiqc_vep(multiqc_file, vep_data_file)
}
vep_data <- read_tsv(vep_data_file)

filter_known <- function(tb, dbsnp) {
  if (is.character(dbsnp)) {
    dbsnp <- read_csv(dbsnp_file)
  }
  tb <- tb |> filter(!is.na(Existing_variation))
  cur_cols <- colnames(tb) |> setdiff("Existing_variation")
  cosmic <- tb |>
    filter(grepl("COSV", Existing_variation)) |>
    select(-Existing_variation)
  tb$Existing_variation <- lapply(tb$Existing_variation, \(x) {
    str_split_1(x, ";")
  })
  in_dbsnp <- tb |>
    filter(map_lgl(Existing_variation, \(x) length(intersect(x, dbsnp$ID)) > 0)) |>
    select(-Existing_variation)
  bind_rows(in_dbsnp, cosmic) |> distinct()
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
