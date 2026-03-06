# [2026-03-06 Fri] you wanna have a function that adds arbitrary subtype labels

library(tidyverse)
library(paletteer)
library(ggpubr)
library(ggplot2)


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
) |>
  distinct() |>
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

main <- ggplot(dummy, aes(x = sample, y = gene, fill = variant)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_paletteer_d("khroma::light")

add_sample_label_plot <- function(
  main_plot,
  label_tb,
  title,
  palette = "ltc::crbhits"
) {
  to_add <- ggplot(label_tb, aes(x = sample, y = 1, fill = label)) +
    geom_tile() +
    scale_fill_paletteer_d(palette) +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      axis.text.y = element_blank()
    ) +
    guides(fill = guide_legend(title = title))
  ggarrange(
    main_plot +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
      ),
    to_add,
    nrow = 2,
    heights = c(10, 1),
    align = "hv"
  )
}

m1 <- add_sample_label_plot(main, subtypes$moffitt, "Moffitt")
## m2 <- add_sample_label_plot(m1, subtypes$kras_subtype, "KRAS subtype")
# TODO: this works, but you gotta add them in a single call

## mof <- ggplot(subtypes$moffitt, aes(x = sample, y = 1, fill = subtype)) +
##   geom_tile() +
##   geom_
