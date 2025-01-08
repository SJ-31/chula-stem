library(tidyverse)
library(ggplot2)
library(glue)
library(here)

here::i_am("./analyses/pdac_vaf_plot.r")
source(here("src", "R", "utils.R"))

vaf_merged_file <- here("analyses", "output", "pdac_vaf_merged.tsv")

if (!file.exists(vaf_merged_file)) {
  data_path <- "/data/project/stemcell/shannc/output/PDAC"
  files <- list.files(data_path, pattern = "8-P[0-9_]+-VEP_small.tsv", recursive = TRUE, full.names = TRUE)
  tsvs <- lapply(files, \(x) {
    tb <- read_with_filename(x, "sample") |> select(sample, VAF, SYMBOL, CLIN_SIG, Consequence)
    tb |>
      group_by(SYMBOL) |>
      summarise(
        VAF = mean(VAF),
        across(where(is.character), \(x) flatten_by(x, collapse = TRUE, unique = FALSE))
      )
  })
  merged <- bind_rows(tsvs) |> mutate(sample = str_extract(sample, "P[0-9]+"))
  write_tsv(merged, vaf_merged_file)
} else {
  merged <- read_tsv(vaf_merged_file) |>
    mutate(across(c(Consequence, CLIN_SIG), into_char_list))
}

## * Format

mode_and_capitalize <- function(char_vec, ignore = c()) {
  if (length(ignore) != 0) {
    char_vec <- discard(char_vec, \(x) x %in% ignore)
  }
  if (length(char_vec) > 1) {
    val <- modes(char_vec) |> first()
  } else {
    val <- first(char_vec)
  }
  val |>
    str_to_title() |>
    str_replace_all("_", " ")
}

min_samples <- 16 # Genes must be found in this number of samples
formatted <- merged |>
  dplyr::filter(
    !is.na(SYMBOL) &
      VAF >= 0.3 &
      !is.na(CLIN_SIG) &
      !is.na(Consequence) &
      map_lgl(Consequence, \(x) "missense_variant" %in% x) &
      map_lgl(CLIN_SIG, \(x) "pathogenic" %in% x)
  ) |>
  group_by(SYMBOL) |>
  filter(n() >= min_samples) |>
  ungroup() |>
  mutate(
    Type = map_chr(Consequence, mode_and_capitalize),
    ClinVar = map_chr(CLIN_SIG, \(x) mode_and_capitalize(x, "benign")),
    ClinVar = case_match(ClinVar, "Na" ~ "Not provided", .default = ClinVar)
  ) |>
  arrange(desc(VAF))

## ** Plot

heatmap_theme <- function() {
  theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    ylab("Gene") +
    xlab("Sample")
}

with_consequence <- formatted |>
  ggplot(aes(x = sample, y = SYMBOL, alpha = VAF, fill = Type)) +
  geom_tile() +
  heatmap_theme()

with_clinsig <- formatted |>
  ggplot(aes(x = sample, y = SYMBOL, alpha = VAF, fill = ClinVar)) +
  geom_tile() +
  heatmap_theme()

blank <- formatted |>
  ggplot(aes(x = sample, y = SYMBOL, alpha = VAF)) +
  geom_tile() +
  heatmap_theme()

ggsave(here("analyses", "output", "pdac_vaf_consequence.png"), with_consequence, dpi = 500)
ggsave(here("analyses", "output", "pdac_vaf_clinsig.png"), with_clinsig, dpi = 500)
ggsave(here("analyses", "output", "pdac_vaf_blank.png"), blank, dpi = 500)
