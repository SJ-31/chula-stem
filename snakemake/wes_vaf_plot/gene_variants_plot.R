library(tidyverse)
library(glue)
library(ggplot2)
library(patchwork)
library(paletteer)

rconfig <- snakemake@config[[snakemake@rule]]
height <- rconfig$height %||% 12
width <- rconfig$width %||% 15
palette_c <- rconfig$palette_c %||% "ggthemes::Green"
palette_d <- rconfig$palette_d %||% "RColorBrewer::Set2"
dpi <- rconfig$dpi
wanted <- snakemake@params$wanted
output <- snakemake@output[[1]]

tb <- read_tsv(
  snakemake@input[[1]]
) |>
  filter(SYMBOL == wanted)

grouped <- tb |>
  group_by(subject, HGVSg) |>
  summarise(
    n_callers = map_dbl(list(unique(INFO_SOURCE)), length),
    mean_af = mean(AF, na.rm = TRUE),
    HGVSp = first(HGVSp),
    POS = first(POS),
    HGVSc = first(HGVSc),
    Consequence = paste0(unique(Consequence), collapse = "&"),
    .groups = "keep"
  ) |>
  mutate(
    HGVSp = ifelse(is.na(HGVSp), NA, str_remove(HGVSp, ".*:")),
    HGVSg = ifelse(is.na(HGVSg), NA, str_remove(HGVSg, ".*:")),
    HGVSc = ifelse(is.na(HGVSc), NA, str_remove(HGVSc, ".*:")),
    var_display = case_when(
      !is.na(HGVSp) ~ HGVSp,
      !is.na(HGVSc) ~ HGVSc,
      .default = HGVSg
    )
  )

ordered_variants <- grouped |>
  ungroup() |>
  select(var_display, POS) |>
  distinct() |>
  deframe() |>
  sort() |>
  names()

var_names2consequence <- grouped |>
  ungroup() |>
  select(var_display, Consequence) |>
  distinct() |>
  deframe() |>
  str_replace_all("_", " ") |>
  str_replace_all("&", ", ") |>
  str_to_title() |>
  str_replace_all("Utr", "UTR") |>
  str_replace_all("3 Prime", "3'") |>
  str_replace_all("5 Prime", "5'")

consequence_axis <- ggplot(
  grouped,
  aes(y = factor(var_display, levels = ordered_variants))
) +
  scale_y_discrete(labels = \(x) var_names2consequence[x]) +
  theme_void() +
  theme(axis.text.y = element_text(hjust = 0))


plot <- ggplot(
  grouped,
  aes(x = subject, y = factor(var_display, levels = ordered_variants))
) +
  geom_tile(
    height = 0.80,
    width = 0.80,
    size = 2,
    aes(fill = mean_af, color = as.factor(n_callers))
  ) +
  theme(
    axis.text.x = element_text(angle = 90),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    legend.title = element_text(face = "bold")
  ) +
  ylab("Variant") +
  guides(
    fill = guide_legend("Mean AF"),
    color = guide_legend("Number of callers")
  ) +
  xlab("Sample") +
  ggtitle(glue("{wanted} variants")) +
  scale_fill_paletteer_c(palette_c) +
  scale_color_paletteer_d(palette_d)

combined <- plot +
  consequence_axis +
  plot_layout(widths = c(10, 0.1), guides = "collect")

ggsave(output, combined, width = width, height = height, dpi = dpi)
