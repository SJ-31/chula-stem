library(tidyverse)
library(here)
suppressWarnings(box::use(
  R / utils[read_existing, alr_normalize],
  dplyr[select, filter]
))

## * GATA6

# Use ALR on the raw counts
housekeeping <- c("MLF2", "SF3B1", "GNB1", "CTBP1", "MYH9")
# Took the top 5 genes in pancreas from HRT atlas, sorting by lowest STD

h2 <- c("FOXH1", "TMEM269", "CD2BP2", "KRT8P33", "RHOQP1", "FREY1")

gata6_hk_plot <- alr_normalize(
  dds[, str_detect(colnames(dds), "PHcase")],
  housekeeping = housekeeping,
  target = "GATA6"
) |>
  ggplot(aes(x = sample, y = housekeeping, fill = n_expr)) +
  geom_tile() +
  geom_text(aes(label = rank), color = "white", size = 6) +
  xlab("Sample") +
  ylab("Housekeeping gene") +
  guides(fill = guide_legend("Normalized expression\n(additive log ratio)")) +
  labs(
    title = "GATA6 Expression"
  )
gata6_hk_plot

ggsave(
  here(workdir, "gata6_hk.png"),
  gata6_hk_plot,
  height = 12,
  width = 12,
  dpi = 500
)
