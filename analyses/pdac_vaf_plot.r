library(tidyverse)
library(ggplot2)
library(glue)
library(here)

here::i_am("./analyses/pdac_vaf_plot.r")
source(here("src", "R", "utils.R"))

vaf_merged_file <- here("analyses", "output", "pdac_vaf_merged.tsv")
vaf_merged_transcripts_file <- here("analyses", "output", "pdac_vaf_merged_transcripts.tsv")

if (!file.exists(vaf_merged_file)) {
  data_path <- "/data/project/stemcell/shannc/output/PDAC"
  files <- list.files(data_path, pattern = "8-P[0-9_]+-VEP_small.tsv$", recursive = TRUE, full.names = TRUE)
  tsvs <- lapply(files, \(x) {
    tb <- read_with_filename(x, "sample") |> select(sample, VAF, Feature, SYMBOL, CLIN_SIG, Consequence)
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
        across(where(is.character), \(x) flatten_by(x, collapse = TRUE, unique = FALSE))
      ) |>
      mutate(census_mutations = map_chr(census_mutations, \(x) {
        if (is.na(x)) {
          return(x)
        }
        first(modes(str_split_1(x, ";") |> discard(\(x) x == "NA")))
      }))
  })
  merged <- bind_rows(tsvs) |> mutate(sample = str_extract(sample, "P[0-9]+"))
  write_tsv(merged, vaf_merged_file)
  q()
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

prettify <- function(tb) {
  tb |> mutate(
    Type = map_chr(Consequence, mode_and_capitalize),
    ClinVar = map_chr(CLIN_SIG, \(x) mode_and_capitalize(x, "benign")),
    ClinVar = case_match(ClinVar, "Na" ~ "Not provided", .default = ClinVar)
  )
}

sum_vafs <- merged |>
  group_by(SYMBOL) |>
  summarise(VAF = sum(VAF)) |>
  tb2map("SYMBOL", "VAF", FALSE) |>
  sort(decreasing = TRUE)

merged$SYMBOL <- factor(merged$SYMBOL, levels = names(sum_vafs))

min_samples <- 16 # Genes must be found in this number of samples


custom_filtered <- merged |>
  dplyr::filter(
    !is.na(SYMBOL) &
      VAF >= 0.3 &
      !is.na(CLIN_SIG) &
      !is.na(Consequence) &
      map_lgl(Consequence, \(x) "missense_variant" %in% x) &
      map_lgl(CLIN_SIG, \(x) "pathogenic" %in% x)
  ) |>
  group_by(SYMBOL) |>
  dplyr::filter(n() >= min_samples) |>
  ungroup()


## * Plot

vaf_heatmap <- function(plot) {
  plot +
    geom_tile() +
    xlab("Sample") + ylab("Gene") +
    theme_grey() +
    theme(
      plot.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.ticks = element_line(size = 0.5),
      axis.text.x = element_text(angle = 90),
    )
}

## ** Custom filters

with_consequence <- custom_filtered |>
  prettify() |>
  ggplot(aes(x = sample, y = SYMBOL, alpha = VAF, fill = Type)) |>
  vaf_heatmap()


with_clinsig <- custom_filtered |>
  prettify() |>
  ggplot(aes(x = sample, y = SYMBOL, alpha = VAF, fill = ClinVar)) |>
  vaf_heatmap()

blank <- custom_filtered |>
  prettify() |>
  ggplot(aes(x = sample, y = SYMBOL, alpha = VAF)) |>
  vaf_heatmap()

ggsave(here("analyses", "output", "pdac_vaf_consequence.png"), with_consequence, dpi = 500)
ggsave(here("analyses", "output", "pdac_vaf_clinsig.png"), with_clinsig, dpi = 500)
ggsave(here("analyses", "output", "pdac_vaf_blank.png"), blank, dpi = 500)

## ** Replicate plot
# Replicate the figure provided
target_genes <- c("KRAS", "TP53", "MUC5B", "KMT2C", "ARID1A", "SMAD4", "GLI3", "CDKN2A")
replicate_figure <- merged |> dplyr::filter(SYMBOL %in% target_genes)

replicate_plot <- replicate_figure |>
  prettify() |>
  ggplot(aes(x = sample, y = SYMBOL, alpha = VAF, fill = Type)) |>
  vaf_heatmap() +
  geom_text(aes(label = round(VAF, 2)), alpha = 1) +
  scale_y_discrete(limits = rev)

ggsave(here("analyses", "output", "pdac_vaf_replicate.png"), replicate_plot, dpi = 500, width = 15)


## ** Cross-referenced
intogen <- read_tsv(here("analyses/data/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"))
census <- read_csv(here("analyses/data/census.csv")) |>
  rename_with(\(x) str_replace_all(x, " ", "_") |> str_remove_all("[()]")) |>
  mutate(Role_in_Cancer = map_chr(Role_in_Cancer, \(x) {
    if (is.na(x)) {
      return(x)
    }
    str_split_1(x, ",") |>
      trimws() |>
      sort() |>
      paste0(collapse = ", ")
  }))

## *** Census

with_census <- inner_join(merged, census, by = join_by(
  x$SYMBOL == y$Gene_Symbol,
  x$census_mutations == y$Mutation_Types
)) |>
  dplyr::filter((Tier == 1) &
    (Hallmark == "Yes") &
    (!is.na(Tumour_TypesSomatic)))


mutation_types <- flatten_by(census$Mutation_Types, ",", collapse = FALSE) |>
  unlist() |>
  map_chr(\(x) trimws(x)) |>
  unique()
# Mis = missense
# T = translocation
# D = large deletion
# F = frameshift
# N = ???
# O = other
# A = amplification
# S = splice site
# M = mesenchymal???

# <2025-01-10 Fri> Genes shown in this plot meet the following criteria
# - have the same consequences as observed in the COSMIC census
#     from Missense, Frameshift and Deletion
# - are Tier 1, so have documented cancer activity
# - Is associated with a known somatic tumor type
# - Has a known hallmark

with_census_plot <- with_census |>
  prettify() |>
  ggplot(aes(
    x = sample, y = factor(SYMBOL, levels = names(sum_vafs)),
    fill = Role_in_Cancer, alpha = VAF
  )) |>
  vaf_heatmap() + geom_text(aes(label = round(VAF, 2))) +
  scale_y_discrete(limits = rev)

ggsave(here("analyses", "output", "pdac_cosmic_census.png"), with_census_plot, dpi = 500, width = 15)


## *** Intogen

intogen_filtered <- intogen |> filter((IS_DRIVER == "TRUE") & (ROLE != "ambiguous") &
  (TOTAL_SAMPLES >= 50))
with_intogen <- merged |>
  separate_longer_delim("Feature", ";") |>
  inner_join(intogen_filtered, by = join_by(SYMBOL, x$Feature == y$TRANSCRIPT)) |>
  group_by(SYMBOL, sample) |>
  summarise(VAF = mean(VAF))
