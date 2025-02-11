library(ggpubr)
source(here::here("analyses", "pdac", "main.R"))

merged <- read_tsv(vaf_merged_file) |>
  mutate(across(c(Consequence, CLIN_SIG), utils$into_char_list))
vep_data <- read_tsv(vep_data_file)


## * Format

mode_and_capitalize <- function(char_vec, ignore = c()) {
  if (length(ignore) != 0) {
    char_vec <- discard(char_vec, \(x) x %in% ignore)
  }
  if (length(char_vec) > 1) {
    val <- utils$modes(char_vec) |> first()
  } else {
    val <- first(char_vec)
  }
  val |>
    str_to_title() |>
    str_replace_all("_", " ")
}

IGNORED_VARIANTS <- c(
  "non_coding_transcript_variant",
  "non_coding_transcript_exon_variant",
  "incomplete_terminal_codon_variant"
)

recode_var_types <- function(types) {
  case_when(
    grepl("utr variant", types) ~ "UTR variant",
    grepl("Splice", types) ~ "Splice variant",
    grepl("Regulatory", types) ~ "Regulatory variant",
    grepl("Inframe", types) ~ "Inframe deletion/insertion",
    grepl("Downstream|Upstream|Intergenic", types) ~ "Intergenic variant",
    grepl("Nmd", types) ~ "NMD transcript variant",
    grepl("mirna|Start|Stop", types) ~ "Other",
    types == "Protein altering variant" ~ "Missense variant",
    types == "Coding sequence variant" ~ "Other",
    .default = types
  )
}

tmb_merged <- vep_data |>
  filter(!key %in% IGNORED_VARIANTS) |>
  filter(key %in% ACCEPTED_CONSEQUENCE) |>
  mutate(
    sample = str_extract(sample, "P[0-9_]+"),
    key = map_chr(key, \(x) str_replace_all(str_to_title(x), "_", " "))
    ## key = recode_var_types(key)
  ) |>
  dplyr::filter(category == "Consequences (all)")

TYPE_ORDER <- tmb_merged |>
  group_by(key) |>
  summarise(mean = mean(value)) |>
  arrange(mean) |>
  pluck("key")

prettify <- function(tb) {
  tb |> mutate(
    type = map_chr(Consequence, \(x) mode_and_capitalize(x, IGNORED_VARIANTS)),
    clinvar = map_chr(CLIN_SIG, \(x) mode_and_capitalize(x, "benign")),
    clinvar = case_match(clinvar, "na" ~ "not provided", .default = clinvar)
    ## type = factor(recode_var_types(type), levels = TYPE_ORDER)
  )
}

merged <- merged |>
  group_by(SYMBOL) |>
  mutate(sum_VAF = sum(VAF), mean_VAF = mean(VAF)) |>
  ungroup()

sum_vafs <- select(merged, SYMBOL, sum_VAF) |>
  distinct() |>
  utils$tb2map("SYMBOL", "sum_VAF", FALSE) |>
  sort(decreasing = TRUE) # Sort order by the summed VAF

by_frequency <- select(merged, SYMBOL, sample) |>
  distinct() |>
  group_by(SYMBOL) |>
  summarise(count = n()) |>
  utils$tb2map("SYMBOL", "count", FALSE) |>
  sort(decreasing = TRUE)

## merged$SYMBOL <- factor(merged$SYMBOL, levels = names(sum_vafs))
merged$SYMBOL <- factor(merged$SYMBOL, levels = names(by_frequency))

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


## * Plot

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

## ** replicate plot
# --- CODE BLOCK ---
sbs <- read_tsv(sbs_merged_file)
# replicate the figure provided
target_genes <- c("KRAS", "TP53", "MUC5B", "KMT2C", "ARID1A", "SMAD4", "GLI3", "CDKN2A")

replicate_figure <- merged |>
  dplyr::filter(SYMBOL %in% target_genes)
n_samples <- sbs$sample |>
  unique() |>
  length()

## *** Add filter
## --- CODE BLOCK ---
filter_version <- "CONSEQUENCE"
if (filter_version == "CLIN_SIG") {
  accepted_clinsig <- c(
    "pathogenic", "likely_pathogenic", "association"
  )
  replicate_figure <- replicate_figure |>
    filter(map_lgl(CLIN_SIG, \(x) length(intersect(x, accepted_clinsig)) > 0))
} else if (filter_version == "CONSEQUENCE") {
  replicate_figure <- filter_known(replicate_figure, dbsnp_file)
  replicate_figure <- mutate(replicate_figure,
    Consequence = lapply(Consequence, \(x) intersect(x, ACCEPTED_CONSEQUENCE))
  ) |>
    filter(map_lgl(Consequence, \(x) length(x) > 0))
} else if (filter_version == "KNOWN") {
  replicate_figure <- filter_known(replicate_figure, dbsnp_file)
  print("Filtering by known variants")
} else if (filter_version == "MUTECT") {
  replicate_figure <- replicate_figure |> filter(grepl("mutect2", SOURCE))
}
order <- replicate_figure$SYMBOL |>
  table() |>
  sort(decreasing = TRUE) |>
  names()
replicate_figure$SYMBOL <- factor(replicate_figure$SYMBOL, levels = order)

## --- CODE BLOCK ---

sample_freq <- replicate_figure |>
  distinct(SYMBOL, sample, .keep_all = TRUE) |>
  group_by(SYMBOL) |>
  summarise(freq = (round(n() / n_samples, 2) * 100) |> as.character() %>% paste0(., " %"))

sbs_plot <- sbs |> ggplot(aes(x = sample, y = count, fill = type)) +
  geom_bar(position = "fill", stat = "identity") +
  guides(fill = guide_legend(title = element_blank())) +
  theme_void() +
  scale_fill_paletteer_d("ggthemes::excel_Depth")

rep_theme <- "tidyquant::tq_light"
axis_title_size <- 12

## *** TMB plot

tmb_max <- group_by(tmb_merged, sample) |>
  summarise(value = sum(value)) |>
  pluck("value") |>
  max()

tmb_plot <- tmb_merged |>
  ggplot(aes(
    x = sample, y = value,
    fill = factor(key, levels = TYPE_ORDER)
  )) +
  geom_bar(position = "stack", stat = "identity") +
  theme_void() +
  ylab("Variant count") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      angle = 90, face = "italic", size = axis_title_size,
      margin = margin(1, 0.5, 0, 0, unit = "cm")
    ),
    axis.ticks.length.y.left = unit(5, "points"),
    axis.line.y = element_line(),
    axis.ticks.y = element_line(), axis.text.y = element_text()
  ) +
  guides(fill = guide_legend(title = element_blank(), position = "bottom")) +
  scale_fill_paletteer_d(rep_theme, drop = FALSE) +
  scale_y_continuous(limits = c(0, tmb_max), expand = c(0, 0), breaks = c(0, tmb_max))

## *** Counts plot
counts_plot <- replicate_figure |>
  distinct(sample, SYMBOL, .keep_all = TRUE) |>
  prettify() |>
  ggplot(aes(y = SYMBOL, fill = factor(type, levels = TYPE_ORDER))) +
  geom_bar() +
  scale_y_discrete(limits = rev) +
  theme_void() +
  theme(
    axis.text.x = element_text(),
    axis.title.x = element_text(face = "italic", size = axis_title_size),
    axis.line.x = element_line(), axis.ticks.x = element_line(),
    axis.ticks.length.x = unit(5, "points")
  ) +
  scale_x_continuous(
    position = "top", limits = c(0, n_samples),
    breaks = c(0, n_samples), expand = c(0, 0)
  ) +
  xlab("Number of samples") +
  guides(fill = "none") +
  scale_fill_paletteer_d(rep_theme, drop = FALSE)

## *** Heatmap
r1 <- replicate_figure |>
  prettify() |>
  ggplot(aes(
    x = sample, y = SYMBOL, # turned off alpha, was confusing
    fill = factor(type, levels = TYPE_ORDER)
  )) |>
  vaf_heatmap() +
  scale_y_discrete(limits = rev) +
  scale_fill_paletteer_d(rep_theme, drop = FALSE) +
  guides(fill = "none") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = axis_title_size),
    axis.title.y = element_text(size = axis_title_size, face = "italic")
  )

r2 <- r1 + scale_y_discrete(limits = rev, labels = rev(sample_freq$freq)) +
  guides(fill = "none") + theme(axis.ticks.y = element_blank())

freq_plot <- ggdraw(get_y_axis(r2))
r3 <- ggdraw(insert_yaxis_grob(r1, get_y_axis(r2),
  position = "right",
  width = unit(0.03, "null")
))

## ** Arranging

# Get and remove legends
type_legend <- ggpubr::get_legend(tmb_plot)
## vaf_legend <- ggpubr::get_legend(r1)
tmb_plot <- tmb_plot + guides(fill = "none")

## legends <- ggarrange(ggdraw(type_legend), NULL, ggdraw(vaf_legend),
##   nrow = 1, widths = c(0.5, -0.3, 0.5)
## )

r1 <- r1 + guides(alpha = "none")

replicate_plot <- ggarrange(
  tmb_plot, NULL, NULL, NULL,
  NULL, NULL, NULL, NULL,
  r1, NULL, freq_plot, counts_plot,
  sbs_plot, NULL, NULL, NULL,
  ncol = 4, nrow = 4, align = "hv",
  widths = c(0.7, 0.01, -0.06, 0.3), heights = c(0.2, -0.08, 0.6, 0.2)
)

final_rep <- ggarrange(replicate_plot, type_legend, ncol = 1, heights = c(0.9, 0.1))
final_rep

save_fn(final_rep, "pdac_vaf_replicate.png")
# --- CODE BLOCK ---
## * Plot metrics
smm <- replicate_figure |>
  group_by(SYMBOL) |>
  summarise(
    median = median(Alt_depth), max = max(Alt_depth), min = min(Alt_depth),
    mean_vaf = median(VAF)
  )

metric_plot <- replicate_figure |>
  group_by(SYMBOL) |>
  mutate(mean_vaf = mean(VAF)) |>
  ggplot(aes(x = SYMBOL, y = log(Alt_depth))) +
  geom_boxplot(aes(fill = mean_vaf)) +
  geom_label(
    data = smm, aes(x = SYMBOL, y = log(median), label = median),
  ) +
  geom_label(data = smm, aes(x = SYMBOL, y = log(min), label = min)) +
  geom_label(data = smm, aes(x = SYMBOL, y = log(max), label = max)) +
  scale_fill_paletteer_c("ggthemes::Green") +
  guides(fill = guide_legend(title = "Mean VAF")) +
  ylab("Log Alternate allele depth") +
  xlab("Gene")
save_fn(metric_plot, "gene_metrics.png")
