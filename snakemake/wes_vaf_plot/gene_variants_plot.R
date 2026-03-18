library(tidyverse)
library(glue)
library(ggplot2)
library(patchwork)
library(ggtext)
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
  filter(SYMBOL == wanted) |>
  mutate(HGVSg = str_remove(HGVSg, "^chr"))

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

format_consequence <- function(vec) {
  str_replace_all(vec, "_", " ") |>
    str_replace_all("&", ", ") |>
    str_to_title() |>
    str_replace_all("Utr", "UTR") |>
    str_replace_all("3 Prime", "3'") |>
    str_replace_all("5 Prime", "5'")
}

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
  format_consequence()

prop_map <- grouped |>
  ungroup() |>
  count(var_display) |>
  mutate(n = {
    mul <- (n / length(unique(grouped$subject))) * 100
    mul |>
      round(2) |>
      as.character() |>
      paste0("%")
  }) |>
  deframe()

color_map <- setNames(
  paletteer_d("Polychrome::dark")[seq_along(unique(var_names2consequence))],
  unique(var_names2consequence)
)

style_color <- function(vec, color) {
  paste0("<b style='color:", color, "'>", vec, "</b>")
}

consequence_axis <- ggplot(
  grouped,
  aes(y = factor(var_display, levels = ordered_variants))
) +
  scale_y_discrete(
    labels = \(x) {
      mapped <- var_names2consequence[x]
      with_color <- style_color(mapped, color_map[mapped])
      paste0(with_color, " (", prop_map[x], ")")
    }
  ) +
  theme_void() +
  theme(axis.text.y = element_markdown(hjust = 0))


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
    axis.ticks.y.right = element_line(),
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
  plot_layout(widths = c(10, 0.5), guides = "collect")

ggsave(output, combined, width = width, height = height, dpi = dpi)
