library(ggpubr)
source(here::here("analyses", "pdac", "main.R"))

merged <- read_tsv(vaf_merged_file) |>
  mutate(across(c(Consequence, CLIN_SIG), utils$into_char_list))
vep_data <- read_tsv(vep_data_file)

var2cons <- read_tsv(vaf_merged_file) |>
  select(Existing_variation, Consequence, CLIN_SIG) |>
  separate_longer_delim(Existing_variation, ";") |>
  distinct(Existing_variation, .keep_all = TRUE)

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
  "incomplete_terminal_codon_variant",
  "intron_variant",
  "3_prime_UTR_variant"
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

prettify <- function(tb, filter_version = "") {
  if (filter_version != "CURATED") {
    to_ignore <- IGNORED_VARIANTS
  } else {
    to_ignore <- c()
  }
  tb |>
    mutate(
      type = map_chr(
        Consequence,
        \(x) mode_and_capitalize(x, to_ignore)
      ),
      clinvar = map_chr(CLIN_SIG, \(x) mode_and_capitalize(x, "benign")),
      clinvar = case_match(clinvar, "na" ~ "not provided", .default = clinvar)
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
    xlab("Sample") +
    ylab("Gene") +
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
# %%

sbs <- read_tsv(sbs_merged_file)
# replicate the figure provided
target_genes <- c(
  "KRAS",
  "TP53",
  "MUC5B",
  "KMT2C",
  "ARID1A",
  "SMAD4",
  "GLI3",
  "CDKN2A"
)

genes_filtered <- merged |>
  dplyr::filter(SYMBOL %in% target_genes)
n_samples <- sbs$sample |>
  unique() |>
  length()

tmp <- read_csv(here(outdir, "select_genes_VEP.csv"))
replicate_figure <- tmp |>
  filter(apply(tmp, 1, \(row) {
    symbol <- row["SYMBOL"]
    to_pass <- CURATED_VARIANTS[[symbol]]
    hgvsp <- row["HGVSp"]
    hgvsg <- row["HGVSg"]
    (hgvsp %in% to_pass) | (hgvsg %in% to_pass)
  })) |>
  group_by(subject, SYMBOL) |>
  summarise(
    Consequence = list(unique(unlist(lapply(
      Consequence,
      \(s) {
        if (!is.na(s)) {
          str_split_1(s, "&")
        } else {
          ""
        }
      }
    )))),
    CLIN_SIG = list(unique(unlist(lapply(
      CLIN_SIG,
      \(s) {
        if (!is.na(s)) {
          str_split_1(s, "&")
        } else {
          ""
        }
      }
    )))),
    Existing_variation = map_chr(
      Existing_variation,
      \(x) str_replace_all(x, "&", ";")
    ),
    VAF = mean(AF),
    Alt_depth = mean(map_dbl(AD, \(x) {
      if (x == ".") {
        NA
      } else {
        splits <- str_split_1(x, ",")
        as.numeric(splits[length(splits)])
      }
    }))
  ) |>
  dplyr::rename(sample = subject) |>
  distinct() |>
  filter(Alt_depth >= 10)

order <- replicate_figure$SYMBOL |>
  table() |>
  sort(decreasing = TRUE) |>
  names()
replicate_figure$SYMBOL <- factor(replicate_figure$SYMBOL, levels = order)

sample_freq <- replicate_figure |>
  distinct(SYMBOL, sample) |>
  group_by(SYMBOL) |>
  summarise(
    freq_raw = (round(n() / n_samples, 2) * 100),
    freq = freq_raw |>
      as.character() %>%
      paste0(., " %")
  ) |>
  arrange(desc(freq_raw))

sbs_plot <- sbs |>
  ggplot(aes(x = sample, y = count, fill = type)) +
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
    x = sample,
    y = value,
    fill = factor(key, levels = TYPE_ORDER)
  )) +
  geom_bar(position = "stack", stat = "identity") +
  theme_void() +
  ylab("Variant count") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      angle = 90,
      face = "italic",
      size = axis_title_size,
      margin = margin(1, 0.5, 0, 0, unit = "cm")
    ),
    axis.ticks.length.y.left = unit(5, "points"),
    axis.line.y = element_line(),
    axis.ticks.y = element_line(),
    axis.text.y = element_text()
  ) +
  guides(fill = guide_legend(title = element_blank(), position = "bottom")) +
  scale_fill_paletteer_d(rep_theme, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, tmb_max),
    expand = c(0, 0),
    breaks = c(0, tmb_max)
  )

## *** Counts plot
counts_plot <- replicate_figure |>
  distinct(sample, SYMBOL, .keep_all = TRUE) |>
  prettify(filter_version) |>
  ggplot(aes(y = SYMBOL, fill = factor(type, levels = TYPE_ORDER))) +
  geom_bar() +
  scale_y_discrete(limits = rev(sample_freq$SYMBOL)) +
  theme_void() +
  theme(
    axis.text.x = element_text(),
    axis.title.x = element_text(face = "italic", size = axis_title_size),
    axis.line.x = element_line(),
    axis.ticks.x = element_line(),
    axis.ticks.length.x = unit(5, "points")
  ) +
  scale_x_continuous(
    position = "top",
    limits = c(0, n_samples),
    breaks = c(0, n_samples),
    expand = c(0, 0)
  ) +
  xlab("Number of samples") +
  guides(fill = "none") +
  scale_fill_paletteer_d(rep_theme, drop = FALSE)

## *** Heatmap
r1 <- replicate_figure |>
  prettify(filter_version) |>
  ggplot(aes(
    x = sample,
    y = SYMBOL, # turned off alpha, was confusing
    fill = factor(type, levels = TYPE_ORDER)
  )) |>
  vaf_heatmap() +
  scale_y_discrete(limits = rev(sample_freq$SYMBOL)) +
  scale_fill_paletteer_d(rep_theme, drop = FALSE) +
  guides(fill = "none") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = axis_title_size),
    axis.title.y = element_text(size = axis_title_size, face = "italic")
  )

r2 <- r1 +
  scale_y_discrete(limits = rev, labels = rev(sample_freq$freq)) +
  guides(fill = "none") +
  theme(axis.ticks.y = element_blank())

freq_plot <- ggdraw(get_y_axis(r2))
r3 <- ggdraw(insert_yaxis_grob(
  r1,
  get_y_axis(r2),
  position = "right",
  width = unit(0.03, "null")
))

## **** for each variant type

heatmap_helper <- function(tb) {
  tb |>
    unnest(cols = c(Consequence)) |>
    prettify(filter_version) |>
    ggplot(aes(
      x = sample,
      y = SYMBOL,
      fill = VAF,
    )) |>
    vaf_heatmap() +
    scale_y_discrete(limits = rev(sample_freq$SYMBOL)) +
    scale_fill_paletteer_c("ggthemes::Blue") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = axis_title_size),
      axis.title.y = element_text(size = axis_title_size, face = "italic")
    ) +
    facet_wrap(~Consequence)
}

with_separate_cons <- heatmap_helper(replicate_figure)
save_fn(with_separate_cons, "pdac_manual_review_separate.png")

## ** Arranging

# Get and remove legends
type_legend <- ggpubr::get_legend(tmb_plot)
## vaf_legend <- ggpubr::get_legend(r1)
tmb_plot <- tmb_plot + guides(fill = "none")

r1 <- r1 + guides(alpha = "none")

replicate_plot <- ggarrange(
  tmb_plot,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  r1,
  NULL,
  freq_plot,
  counts_plot,
  sbs_plot,
  NULL,
  NULL,
  NULL,
  ncol = 4,
  nrow = 4,
  align = "hv",
  widths = c(0.7, 0.01, -0.06, 0.3),
  heights = c(0.2, -0.08, 0.6, 0.2)
)

final_rep <- ggarrange(
  replicate_plot,
  type_legend,
  ncol = 1,
  heights = c(0.9, 0.1)
)
final_rep

save_fn(final_rep, "pdac_vaf_replicate.png")
# --- CODE BLOCK ---
## * Plot metrics
smm <- replicate_figure |>
  group_by(SYMBOL) |>
  summarise(
    median = median(Alt_depth),
    max = max(Alt_depth),
    min = min(Alt_depth),
    mean_vaf = median(VAF)
  )

metric_plot <- replicate_figure |>
  group_by(SYMBOL) |>
  mutate(mean_vaf = mean(VAF)) |>
  ggplot(aes(x = SYMBOL, y = log(Alt_depth))) +
  geom_boxplot(aes(fill = mean_vaf)) +
  geom_label(
    data = smm,
    aes(x = SYMBOL, y = log(median), label = median),
  ) +
  geom_label(data = smm, aes(x = SYMBOL, y = log(min), label = min)) +
  geom_label(data = smm, aes(x = SYMBOL, y = log(max), label = max)) +
  scale_fill_paletteer_c("ggthemes::Green") +
  guides(fill = guide_legend(title = "Mean VAF")) +
  ylab("Log Alternate allele depth") +
  xlab("Gene")
save_fn(metric_plot, "gene_metrics.png")

## * Summary of known variants

source(here("src", "R", "utils.R"))

clinsig_hierarchy <- c(
  "not_provided" = 0,
  "pathogenic" = 3,
  "association" = 1,
  "likely_pathogenic" = 2
)

get_existing_var <- function(tb, symbols) {
  filter(tb, SYMBOL %in% symbols) |>
    separate_longer_delim(Existing_variation, ";") |>
    select(-Consequence, -CLIN_SIG) |>
    inner_join(var2cons) |>
    mutate(
      clinsig = map_chr(CLIN_SIG, \(chars) {
        if (!is.na(chars)) {
          splits <- str_split_1(chars, ";")
          sorted <- clinsig_hierarchy[splits] |> sort(decreasing = TRUE)
          names(sorted[1])
        } else {
          ""
        }
      }),
      Consequence = map_chr(
        Consequence,
        \(vec) flatten_by(vec, unique = TRUE, collapse = TRUE)
      )
    )
}

targets <- get_existing_var(replicate_figure, target_genes) |>
  select(-where(is.numeric), -CLIN_SIG) |>
  mutate(SYMBOL = as.character(SYMBOL))

## lapply(unique(targets$SYMBOL), \(x) {
##   current <- filter(targets, SYMBOL == x)
##   table(current$Existing_variation) |>
##     as.data.frame() |>
##     write_tsv(here(outdir, glue("{x}_existing_vars.tsv")))
## })

# bcftoo
# TODO: lookup these existing variants

## * Variant tables

table_outdir <- here(outdir, "vtables")
to_tables <- tmp |>
  filter(apply(tmp, 1, \(row) {
    symbol <- row["SYMBOL"]
    to_pass <- CURATED_VARIANTS[[symbol]]
    hgvsp <- row["HGVSp"]
    hgvsg <- row["HGVSg"]
    (hgvsp %in% to_pass) | (hgvsg %in% to_pass)
  })) |>
  mutate(id = case_when(is.na(HGVSp) ~ HGVSg, .default = HGVSp)) |>
  select(subject, id, SYMBOL, AF, AD) |>
  separate_wider_delim(
    cols = AD,
    delim = ",",
    names = c("Ref_depth", "Alt_depth")
  ) |>
  group_by(subject, SYMBOL, id) |>
  summarize(
    id = first(id),
    AF = mean(AF),
    Ref_depth = mean(as.numeric(Ref_depth)),
    Alt_depth = mean(as.numeric(Alt_depth))
  ) |>
  mutate(value = paste0(Ref_depth, ",", Alt_depth, " (", round(AF, 2), ")"))

for (sym in unique(to_tables$SYMBOL)) {
  to_tables |>
    filter(SYMBOL == sym) |>
    select(-SYMBOL) |>
    pivot_wider(
      names_from = id,
      id_cols = subject,
      values_from = value,
      values_fill = "-"
    ) |>
    mutate(across(where(is.double), \(x) round(x, 2))) |>
    write_tsv(here(table_outdir, glue("{sym}.tsv")))
}

## * Caller inconsistencies

tmp |>
  mutate(AD2 = AD) |>
  separate_wider_delim(
    cols = AD2,
    delim = ",",
    names = c("Ref_depth", "Alt_depth"),
    too_few = "align_end",
    too_many = "drop"
  ) |>
  mutate(
    Ref_depth = as.numeric(Ref_depth),
    Alt_depth = as.numeric(Alt_depth)
  ) |>
  filter(AF != 1 & Ref_depth == 0) |>
  select(
    subject,
    INFO_SOURCE,
    SYMBOL,
    HGVSg,
    HGVSp,
    Consequence,
    GT,
    AD,
    AF,
    Existing_variation
  ) |>
  write_tsv(here(outdir, "inconsistent_gt_ad.tsv"))

## * Oncoprint

library(ComplexHeatmap)

oncoprint_allowed <- c(
  "missense_variant" = "blue",
  "frameshift_variant" = "red",
  "downstream_gene_variant" = "green",
  "upstream_gene_variant" = "brown",
  "stop_gained" = "yellow",
  "splice_region_variant" = "#8839ef",
  "inframe_deletion" = "#e64553"
)

to_oncoprint <- replicate_figure |>
  select(sample, SYMBOL, Consequence) |>
  mutate(
    Consequence = lapply(
      Consequence,
      \(x) intersect(x, names(oncoprint_allowed))
    )
  ) |>
  group_by(SYMBOL, sample) |>
  pivot_wider(
    names_from = sample,
    values_from = Consequence,
    values_fn = \(x) paste0(unique(x[[1]]), collapse = ";")
  ) |>
  mutate(across(where(is.character), \(x) replace(x, is.na(x), ""))) |>
  column_to_rownames(var = "SYMBOL") |>
  as.matrix()


offset <- 1 / length(oncoprint_allowed) / length(oncoprint_allowed)
alter_fun <- sapply(
  names(oncoprint_allowed),
  \(var) {
    \(x, y, w, h) {
      grid.rect(
        x = x,
        y = y,
        width = w,
        height = h / length(oncoprint_allowed),
        gp = gpar(fill = oncoprint_allowed[var])
      )
    }
  },
  simplify = FALSE,
  USE.NAMES = TRUE
)

## oncoPrint(
##   to_oncoprint,
##   get_type = \(x) str_split_1(x, ";"),
##   alter_fun = alter_fun,
##   show_column_names = TRUE,
##   col = oncoprint_allowed
## )

## --- CODE BLOCK ---
