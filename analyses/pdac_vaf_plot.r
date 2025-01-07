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

## * Plotting
min_samples <- 17 # Genes must be found in this number of samples
formatted <- merged |>
  filter(
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
  mutate(Type = map_chr(Consequence, \(x) {
    if (length(x) > 1) {
      val <- modes(x) |> first()
    } else {
      val <- first(x)
    }
    val |>
      str_to_title() |>
      str_replace_all("_", " ")
  }))

plot <- formatted |> ggplot(aes(x = sample, y = SYMBOL, alpha = VAF, fill = Type)) +
  geom_tile() +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  ylab("Gene") +
  xlab("Sample")

plot
