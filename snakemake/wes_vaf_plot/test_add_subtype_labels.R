# [2026-03-06 Fri] you wanna have a function that adds arbitrary subtype labels

library(tidyverse)
library(paletteer)
library(checkmate)
library(here)
library(patchwork)
library(ggplot2)

Sys.setenv(R_SRC = here("src", "R"))
source(here("src", "R", "plotting.R"))

var_types <- c(
  "missense",
  "nonsense",
  "frameshift",
  "stop_gaind",
  "deletion",
  NA
)
genes <- c("P53", "KRAS", "MLH1", "APC", "AIM2", "GATA6")
n_sample <- 10
samples <- paste0("P", seq(n_sample))

moffitt <- c("basal-like", "classical")
kras_subtype <- c("wt", "balanced", "minor", "major")

n <- 1000
dummy <- tibble(
  sample = sample(samples, size = n, replace = TRUE),
  gene = sample(genes, size = n, replace = TRUE),
  variant = sample(var_types, size = n, replace = TRUE)
)

var_freq_plot <- ggplot(dummy, aes(x = sample, fill = variant)) +
  geom_bar() +
  theme_void() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  scale_fill_paletteer_d("khroma::light") +
  guides(fill = "none")

dummy_distinct <- dummy |>
  distinct(sample, gene, .keep_all = TRUE) |>
  filter(!is.na(variant))


subtypes <- list(
  moffitt = tibble(
    sample = samples,
    label = sample(moffitt, size = n_sample, replace = TRUE)
  ),
  kras_subtype = tibble(
    sample = samples,
    label = sample(kras_subtype, size = n_sample, replace = TRUE)
  )
)

TILE_CALL <- geom_tile(width = 0.95, height = 0.95)

main <- ggplot(dummy_distinct, aes(x = sample, y = gene, fill = variant)) +
  TILE_CALL +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  scale_fill_paletteer_d("khroma::light")


sample_freq <- dummy_distinct |>
  group_by(gene) |>
  summarise(freq = n() / length(unique(dummy_distinct$sample)))

var_count_plot <- ggplot(dummy_distinct, aes(y = gene, fill = variant)) +
  geom_bar() +
  scale_fill_paletteer_d("khroma::light") +
  scale_y_discrete(labels = \(lims) {
    deframe(sample_freq)[lims]
  }) +
  theme_void() +
  theme(axis.text.y = element_text()) +
  guides(fill = "none")

legend <- ggpubr::get_legend(main)
main <- main + guides(fill = "none")

m1 <- combine_sample_label_plots(
  subtypes,
  palettes = list(kras_subtype = "ltc::crbhits", moffitt = "MetBrewer::Cross")
)

blank_margin <- theme(plot.margin = margin(0, 0, 0, 0, "cm"))

y_plot <- cowplot::ggdraw(cowplot::get_y_axis(var_count_plot))
plots <- list(
  var_freq_plot,
  main,
  m1,
  var_count_plot + theme(axis.text.y = element_blank()),
  legend,
  y_plot + blank_margin
)


## r2 <- main +

## TODO: want to check getting the proportion axis
final <- wrap_plots(plots) +
  plot_layout(
    design = "
1##
264
3##
55#
",
    heights = c(3, 10, 1),
    widths = c(1, 0, 1)
  )
final


ggsave(here("test_subtypes.pdf"), final)

# TODO: this works, but you gotta add them in a single call

## mof <- ggplot(subtypes$moffitt, aes(x = sample, y = 1, fill = subtype)) +
##   geom_tile() +
##   geom_
