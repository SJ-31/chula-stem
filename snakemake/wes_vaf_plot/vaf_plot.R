suppressMessages({
  library(patchwork)
  library(checkmate)
  library(glue)
  library(ggplot2)
  library(tidyverse)
  library(paletteer)
  library(ggtext)
  library(logger)
  library(here)
  utils <- new.env()
  source(here("src", "R", "utils.R"), local = utils)
  source(here("src", "R", "plotting.R"))
  tmb_merged <- read_tsv(snakemake@input$tmb)
  sbs <- read_tsv(snakemake@input$sbs)
  target_genes <- snakemake@params$wanted_genes
  combined_vep <- read_tsv(snakemake@input$combined_vep)
  config <- snakemake@config
  label_spec <- config$sample_labels
  samples <- snakemake@params$samples |> sort()
  nd_label <- config$no_data_label %||% "ND"
  samples_with_wes <- samples
})

print(glue("Available samples: {paste0(samples, collapse = ',')}"))

with_no_data <- c()
TILE_CALL <- geom_tile(width = 0.95, height = 0.95)

# [2026-03-16 Mon] BUG: the samples prefixed with N aren't showing up
# Might be because they have TOOL_SOURCE as the source tag and not "SOURCE". So
# gotta rename them
# More importantly, they use a different system for chr names... Need to relabel chrs
# Add a how_to for both of these

## * Get extra sample labels

if (!is.null(label_spec)) {
  assert_list(label_spec, names = "unique")
  extra_labels <- lapply(label_spec, \(spec) {
    assert_list(spec)
    assert_names(names(spec), must.include = c("palette", "file"))
    tb <- read_tsv(spec$file) |> mutate(label = as.character(label))
    others <- samples |> discard(\(s) s %in% tb$sample)
    new_samples <- tb |>
      pluck("sample") |>
      discard(\(s) s %in% samples)
    if (config$add_label_samples %||% TRUE) {
      samples <<- unique(c(samples, new_samples))
      with_no_data <<- unique(c(with_no_data, new_samples))
    }
    tb |> bind_rows(tb, tibble(sample = others, label = NA))
  })
  label_palettes <- lapply(label_spec, \(s) s$palette)
  extra_labels[["Exome data available"]] <- tibble(sample = samples) |>
    mutate(label = case_when(sample %in% with_no_data ~ "N", .default = "Y"))
  label_palettes[["Exome data available"]] <- "ggsci::alternating_igv"
  samples <<- sort(samples)
  label_plot <- combine_sample_label_plots(
    extra_labels,
    palettes = label_palettes
  ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
    scale_x_discrete(limits = samples)
} else {
  label_plot <- NULL
}


accepted_consequence <- config$accepted_consequence
min_alt_depth <- config$variant_calling$min_alt_depth
CURATED_VARIANTS <- config$allowed_variants
ONLY_CURATED <- config$only_curated

save_plot_smk_params <- function(plot, filename) {
  ggsave(
    filename,
    plot = plot,
    dpi = config$plot$dpi,
    width = config$plot$width,
    height = config$plot$height
  )
}

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


prettify <- function(tb) {
  tb |>
    mutate(
      type = map_chr(
        Consequence,
        \(x) mode_and_capitalize(x)
      ),
      clinvar = map_chr(CLIN_SIG, \(x) mode_and_capitalize(x, "benign")),
      clinvar = replace_values(clinvar, "na" ~ "not provided")
    )
}

## * Plot

vaf_heatmap <- function(plot) {
  plot +
    TILE_CALL +
    xlab("Sample") +
    ylab("Gene") +
    theme_grey() +
    theme(
      plot.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_line(linewidth = 0.5),
      axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
    )
}

n_samples <- samples |>
  unique() |>
  length()

add_groupings <- FALSE

extra <- c()
if ("cases" %in% names(config$extra) && length(config$extra$cases) != 0) {
  add_groupings <- TRUE
  extra <- config$extra$cases
}

if (config$variant_calling$protein_only %||% FALSE) {
  log_info("Filtering by coding variants...")
  log_info("Count before: {nrow(combined_vep)}")
  combined_vep <- filter(combined_vep, !is.na(HGVSp))
  log_info("Count after: {nrow(combined_vep)}")
}

replicate_figure <- combined_vep |>
  filter(apply(combined_vep, 1, \(row) {
    if (ONLY_CURATED) {
      symbol <- row["SYMBOL"]
      to_pass <- CURATED_VARIANTS[[symbol]]
      hgvsp <- row["HGVSp"]
      hgvsg <- row["HGVSg"]
      (hgvsp %in% to_pass) | (hgvsg %in% to_pass)
    } else {
      TRUE
    }
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
    Existing_variation = list(unique(unlist(lapply(
      Existing_variation,
      \(x) str_replace_all(x, "&", ";")
    )))),
    VAF = mean(AF),
    Alt_depth = mean(
      map_dbl(AD, \(x) {
        if (x == ".") {
          NA
        } else {
          splits <- str_split_1(x, ",")
          as.numeric(splits[length(splits)])
        }
      }),
      na.rm = TRUE
    ),
    HGVSp = paste0(HGVSp, collapse = ";")
  ) |>
  dplyr::rename(sample = subject) |>
  distinct()


if (!is.null(min_alt_depth)) {
  replicate_figure <- filter(replicate_figure, Alt_depth >= min_alt_depth)
}


if (!ONLY_CURATED) {
  replicate_figure <- mutate(
    replicate_figure,
    Consequence = lapply(
      Consequence,
      \(csqs) keep(csqs, \(x) x %in% accepted_consequence)
    )
  ) |>
    filter(unlist(lapply(Consequence, length)) >= 1)
  tmb_merged <- filter(tmb_merged, key %in% accepted_consequence)
}

TYPE_ORDER <- tmb_merged |>
  group_by(key) |>
  summarise(mean = mean(value)) |>
  arrange(mean) |>
  pluck("key")
TYPE_ORDER_TITLE <- str_to_title(TYPE_ORDER) |> str_replace_all("_", " ")

order <- replicate_figure$SYMBOL |>
  table() |>
  sort(decreasing = TRUE) |>
  names()
replicate_figure$SYMBOL <- factor(replicate_figure$SYMBOL, levels = order)

sample_freq <- replicate_figure |>
  filter(sample %in% samples_with_wes) |>
  distinct(SYMBOL, sample) |>
  group_by(SYMBOL) |>
  summarise(
    freq_raw = (round(n() / length(samples_with_wes), 2) * 100),
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
  scale_fill_paletteer_d("ggthemes::excel_Depth") +
  scale_x_discrete(limits = samples)

rep_theme <- "tidyquant::tq_light"
axis_title_size <- 12

# Add any extra samples from the labels
if (length(with_no_data) > 0) {
  no_data_tb <- tibble(
    sample = with_no_data,
    Consequence = list(nd_label),
    CLIN_SIG = list(nd_label)
  ) |>
    mutate(
      SYMBOL = lapply(sample, \(.) {
        unique(replicate_figure$SYMBOL)
      })
    ) |>
    unnest(SYMBOL)
  replicate_figure <- bind_rows(replicate_figure, no_data_tb)
}


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
    axis.text.y = element_text(),
    plot.margin = unit(c(30, 0, 0, 0), "points")
  ) +
  guides(fill = guide_legend(title = element_blank(), position = "bottom")) +
  scale_fill_paletteer_d(rep_theme, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, tmb_max),
    expand = c(0, 0),
    breaks = c(0, tmb_max)
  ) +
  scale_x_discrete(limits = samples)


## *** Counts plot
counts_plot <- replicate_figure |>
  ## distinct(sample, SYMBOL, .keep_all = TRUE) |> # [2026-03-09 Mon] Pretty sure you don't need this
  filter(sample %in% samples_with_wes) |>
  prettify() |>
  ggplot(aes(y = SYMBOL, fill = factor(type, levels = TYPE_ORDER_TITLE))) +
  geom_bar() +
  scale_y_discrete(limits = rev(sample_freq$SYMBOL), labels = \(old) {
    paste0(deframe(select(sample_freq, SYMBOL, freq))[old], " ")
  }) +
  theme_void() +
  theme(
    axis.text.x = element_text(),
    axis.title.x = element_text(face = "italic", size = axis_title_size),
    axis.text.y = element_text(),
    axis.line.x = element_line(),
    axis.ticks.x = element_line(),
    axis.ticks.length.x = unit(5, "points"),
    plot.margin = unit(c(0, 30, 0, 0), "points")
  ) +
  scale_x_continuous(
    ## position = "top",
    limits = c(0, length(samples_with_wes)),
    breaks = c(0, length(samples_with_wes)),
    expand = c(0, 0)
  ) +
  xlab("Number of samples (with WES)") +
  guides(fill = "none") +
  scale_fill_paletteer_d(rep_theme, drop = FALSE)

## *** Save plotting data

replicate_figure |>
  unnest(cols = c(Consequence)) |>
  separate_longer_delim(HGVSp, ";") |>
  prettify() |>
  distinct(sample, SYMBOL, type, HGVSp, .keep_all = TRUE) |>
  select(-c(Consequence, CLIN_SIG)) |>
  write_tsv(snakemake@output$plot_data)

## *** Heatmap

r1 <- replicate_figure |>
  prettify() |>
  ggplot(aes(
    x = sample,
    y = SYMBOL,
    fill = factor(type, levels = TYPE_ORDER_TITLE)
  )) |>
  vaf_heatmap() +
  scale_y_discrete(limits = rev(sample_freq$SYMBOL)) +
  scale_fill_paletteer_d(rep_theme, drop = FALSE) +
  guides(fill = "none") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = axis_title_size),
    axis.title.y = element_text(size = axis_title_size, face = "italic")
  ) +
  scale_x_discrete(limits = samples)

if (add_groupings) {
  replicate_figure <- local({
    is_paired <- config$variant_calling$paired
    replicate_figure |>
      mutate(
        group = case_when(
          is_paired & sample %in% extra ~ "tumor-only",
          is_paired & !sample %in% extra ~ "paired",
          !is_paired & sample %in% extra ~ "tumor-only",
          .default = "paired"
        )
      )
  })
  tmp <- replicate_figure |>
    ungroup() |>
    ungroup() |>
    arrange(group, sample) %>%
    mutate(
      labs = glue_data(
        .,
        "<span style='color: {if_else(group == 'tumor-only', 'purple', 'black')}'>{sample}</span>"
      )
    ) |>
    distinct(sample, .keep_all = TRUE)
  labs <- tmp$labs |> `names<-`(tmp$sample)
  sample_names <- tmp$sample
  r1 <- r1 +
    theme(axis.text.x = element_markdown()) +
    scale_x_discrete(labels = labs, limits = samples)
}

## **** for each variant type

heatmap_helper <- function(tb) {
  tb |>
    unnest(cols = c(Consequence)) |>
    prettify() |>
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

with_separate_cons <- heatmap_helper(replicate_figure) +
  scale_x_discrete(limits = samples)
if (add_groupings) {
  with_separate_cons <- with_separate_cons +
    theme(axis.text.x = element_markdown()) +
    scale_x_discrete(labels = labs, limits = samples)
}

## ** Arranging

# Get and remove legends
type_legend <- ggpubr::get_legend(tmb_plot)
## vaf_legend <- ggpubr::get_legend(r1)
tmb_plot <- tmb_plot + guides(fill = "none")

r1 <- r1 + guides(alpha = "none")

# TODO: [2026-03-09 Mon] gotta check this works
# TODO: revert to grabbing the y axis from counts plot if the margins don't work well
to_arrange <- list(
  tmb_plot, # 1
  r1, # 2
  counts_plot, # 3
  sbs_plot, # 4
  type_legend # 5
)
if (!is.null(label_plot)) {
  to_arrange[[2]] <- r1 +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  guide_area <- guide_area()
  label_plot <- label_plot + xlab("Sample")
  to_arrange <- c(to_arrange, label_plot, guide_area)
  design <- "
#1#
723
76#
74#
75#
"
  heights <- c(3, 12, 3, 3, 5)
  widths <- c(1, 10, 3)
} else {
  design <- "
1#
23
35
4#
"
  widths <- c(10, 4)
  heights <- c(3, 10, 3, 5)
}

final_rep <- wrap_plots(
  to_arrange,
  design = design,
  heights = heights,
  widths = widths
) +
  plot_layout(guides = "collect")


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


if (exists("snakemake")) {
  save_plot_smk_params(with_separate_cons, snakemake@output$vaf_separate)
  save_plot_smk_params(final_rep, snakemake@output$vaf)
  save_plot_smk_params(metric_plot, snakemake@output$gene_metrics)
}
