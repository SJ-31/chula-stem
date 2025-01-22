library(tidyverse)
library(ggplot2)
library(glue)
library(here)
library(reticulate)
library(ggpubr)
library(cowplot)

here::i_am("./analyses/pdac/vaf_plot.r")
utils <- new.env()
source(here("src", "R", "utils.R"), local = utils)
py_utils <- new.env()
reticulate::source_python(here("src", "chula_stem", "utils.py"), py_utils)

vaf_merged_file <- here("analyses", "output", "pdac_vaf_merged.tsv")
sbs_merged_file <- here("analyses", "output", "pdac_sbs.tsv")
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
        across(where(is.character), \(x) flatten_by(x, collapse = TRUE, unique = FALSE))
      ) |>
      mutate(census_mutations = map_chr(census_mutations, \(x) {
        if (is.na(x)) {
          return(x)
        }
        first(modes(str_split_1(x, ";") |> discard(\(x) x == "NA")))
      }))
  })
  merged <- bind_rows(tsvs) |> mutate(sample = str_extract(sample, "P[0-9_]+"))
  write_tsv(merged, vaf_merged_file)
  q()
} else {
  merged <- read_tsv(vaf_merged_file) |>
    mutate(across(c(Consequence, CLIN_SIG), into_char_list))
}

multiqc_file <- here(data_path, "8-cohort-MultiQC_data", "vep.txt")
vep_data_file <- here("analyses", "output", "pdac_vep_data.tsv")
if (!file.exists(vep_data_file)) {
  py_utils$parse_multiqc_vep(multiqc_file, vep_data_file)
}
vep_data <- read_tsv(vep_data_file)

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

recode_var_types <- function(types) {
  case_when(
    grepl("utr variant", types) ~ "UTR variant",
    grepl("Stop", types) ~ "Stop variant",
    grepl("Start", types) ~ "Start variant",
    grepl("Splice", types) ~ "Splice variant",
    grepl("Regulatory", types) ~ "Regulatory variant",
    grepl("Inframe", types) ~ "Inframe variant",
    grepl("Non coding", types) ~ "Non-coding sequence variant",
    types == "Protein altering variant" ~ "Missense variant",
    .default = types
  )
}

prettify <- function(tb) {
  tb |> mutate(
    type = map_chr(Consequence, mode_and_capitalize),
    clinvar = map_chr(CLIN_SIG, \(x) mode_and_capitalize(x, "benign")),
    clinvar = case_match(clinvar, "na" ~ "not provided", .default = clinvar),
    type = recode_var_types(type)
  )
}

merged <- merged |>
  group_by(SYMBOL) |>
  mutate(sum_VAF = sum(VAF), mean_VAF = mean(VAF)) |>
  ungroup()

sum_vafs <- select(merged, SYMBOL, sum_VAF) |>
  distinct() |>
  tb2map("SYMBOL", "sum_VAF", FALSE) |>
  sort(decreasing = TRUE)

mean_vafs <- select(merged, mean_VAF, sum_VAF) |>
  distinct() |>
  tb2map("mean_VAF", "sum_VAF", FALSE) |>
  sort(decreasing = TRUE)

# todo: <2025-01-21 tue> this isn't sorting correctly
vaf_levels <- names(mean_vafs) |>
  as.numeric() |>
  round(2) |>
  as.character() |>
  unique()

merged$SYMBOL <- factor(merged$SYMBOL, levels = names(sum_vafs))
merged$mean_vaf_char <- factor(as.character(round(merged$mean_VAF, 2)),
  levels = vaf_levels
)

min_samples <- 16 # genes must be found in this number of samples

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


## * plot

vaf_heatmap <- function(plot) {
  plot +
    geom_tile() +
    xlab("sample") + ylab("gene") +
    theme_grey() +
    theme(
      plot.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.ticks = element_line(size = 0.5),
      axis.text.x = element_text(angle = 90),
    )
}

## ** custom filters

with_consequence <- custom_filtered |>
  prettify() |>
  ggplot(aes(x = sample, y = SYMBOL, alpha = VAF, fill = type)) |>
  vaf_heatmap()


with_clinsig <- custom_filtered |>
  prettify() |>
  ggplot(aes(x = sample, y = SYMBOL, alpha = VAF, fill = clinvar)) |>
  vaf_heatmap()

blank <- custom_filtered |>
  prettify() |>
  ggplot(aes(x = sample, y = SYMBOL, alpha = VAF)) |>
  vaf_heatmap()

ggsave(here("analyses", "output", "pdac_vaf_consequence.png"), with_consequence, dpi = 500)
ggsave(here("analyses", "output", "pdac_vaf_clinsig.png"), with_clinsig, dpi = 500)
ggsave(here("analyses", "output", "pdac_vaf_blank.png"), blank, dpi = 500)


## ** replicate plot
sbs <- read_tsv(sbs_merged_file)
# replicate the figure provided
target_genes <- c("KRAS", "TP53", "MUC5B", "KMT2C", "ARID1A", "SMAD4", "GLI3", "CDKN2A")
replicate_figure <- merged |> dplyr::filter(SYMBOL %in% target_genes)
n_samples <- sbs$sample |>
  unique() |>
  length()

sbs_plot <- sbs |> ggplot(aes(x = sample, y = count, fill = type)) +
  geom_bar(position = "fill", stat = "identity") +
  guides(fill = guide_legend(title = element_blank())) +
  guides(fill = guide_legend(title = element_blank())) +
  theme_void()

# TODO: you don't actually know yet how to get tmb so this is some dummy data
tmb_merged <- vep_data |>
  mutate(
    sample = str_extract(sample, "P[0-9_]+"),
    key = map_chr(key, \(x) str_replace_all(str_to_title(x), "_", " ")),
    key = recode_var_types(key)
  ) |>
  dplyr::filter(category == "Consequences (all)")
tmb_plot <- tmb_merged |> ggplot(aes(x = sample, y = value, fill = key)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_void() +
  theme(axis.title.x = element_blank())

counts_plot <- replicate_figure |>
  prettify() |>
  ggplot(aes(y = SYMBOL, fill = type)) +
  geom_bar() +
  scale_y_discrete(limits = rev) +
  theme_void() +
  theme(
    axis.text.x = element_text(), axis.title.x = element_text(),
    axis.line.x = element_line(), axis.ticks.x = element_line(),
    axis.ticks.length.x = unit(5, "points")
  ) +
  scale_x_continuous(
    position = "top", limits = c(0, n_samples),
    breaks = c(0, n_samples), expand = c(0, 0)
  ) +
  xlab("Number of samples") +
  guides(fill = "none")

r1 <- replicate_figure |>
  prettify() |>
  ggplot(aes(x = sample, y = SYMBOL, alpha = VAF, fill = type)) |>
  vaf_heatmap() +
  scale_y_discrete(limits = rev) + guides(fill = "none")

r2 <- replicate_figure |>
  prettify() |>
  ggplot(aes(x = sample, y = mean_vaf_char, alpha = VAF, fill = type)) |>
  vaf_heatmap() +
  scale_y_discrete(limits = rev)

# Get and remove legends
type_legend <- ggpubr::get_legend(tmb_plot)
vaf_legend <- ggpubr::get_legend(r1)
tmb_plot <- tmb_plot + guides(fill = "none")
r1 <- r1 + guides(alpha = "none")

replicate_plot <- ggdraw(insert_yaxis_grob(r1, get_y_axis(r2, position = "right")))

# <2025-01-22 Wed> Sad hack to change widths
replicate_plot <- ggarrange(
  tmb_plot, NULL, NULL,
  NULL, NULL, NULL,
  r1, NULL, counts_plot,
  sbs_plot, NULL, NULL,
  ncol = 3, nrow = 4, align = "hv",
  widths = c(0.7, -0.07, 0.3), heights = c(0.2, -0.08, 0.6, 0.2)
)
replicate_plot

# TODO: need graph of TMB
# TODO: need sideways bar plot showing sample count
# TODO: need column plot showing the frequency of mutation types

ggsave(here("analyses", "output", "pdac_vaf_replicate.png"),
  replicate_plot,
  dpi = 500, width = 15
)


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
