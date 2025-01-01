library(tidyverse)
library(GenomicRanges)
library(grid)
library(ggplot2)
cns_file <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.call.cns"
# Contains copy number calls
cnr_file <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.cnr"
# Contains the log2 segmentation data
cnn_file <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.targetcoverage.cnn"

facets_rds <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/5-patient_10-Facets/5-patient_10-Facets_hisens.rds"

format_chr_data <- function(tb, divide_by = 1000) {
  chr_levels <- c(as.character(1:22), "X", "Y")
  tb |> mutate(
    start = start / divide_by,
    end = end / divide_by,
    mid = (end - start) / 2 + start,
    chromosome = factor(tb$chromosome, levels = chr_levels)
  )
}

cns <- read_tsv(cns_file) |> format_chr_data()
cnn <- read_tsv(cnn_file)
cnr <- read_tsv(cnr_file) |> format_chr_data()
facets <- readRDS(facets_rds)
antitarget <- cnr |> filter(gene == "Antitarget")
# Off-target bins typically going to be much wider than target bins


tb2map <- function(tb, keys, values) {
  as.list(tb[[values]]) |> `names<-`(tb[[keys]])
}

# <2025-01-01 Wed> Should just represent the log2 values as the midpoint
# of each bin range, as that is what the cnvkit algorithm seems to do anyway
# And should probably only plot target regions
target <- cnr |> filter(gene != "Antitarget")

COLORS <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = TRUE)]
cn_colors <- c("blue", "red", "green", "orange")

color_legend <- legendGrob(as.character(seq(cn_colors) - 1),
  pch = 19,
  gp = gpar(col = cn_colors)
)

chr <- 19


plot_cnvkit <- function(ratios, cn_calls, chr) {
  if (!missing(chr) && length(chr) >= 1) {
    cn_calls <- cn_calls |> filter(chromosome %in% chr)
    ratios <- ratios |> filter(chromosome %in% chr)
  }
  plot <- ratios |>
    inner_join(cn_calls,
      by = join_by(chromosome, between(mid, y$start, y$end)),
      suffix = c("", "_right")
    ) |>
    ggplot(aes(x = mid, y = log2, alpha = weight), stat = "identity") +
    xlab("Position (mb)") +
    ylab("Copy Number Ratio (log2)") +
    geom_point() +
    geom_segment(aes(
      x = start_right, xend = end_right,
      y = log2_right, yend = log2_right,
      color = factor(cn)
    ), alpha = 0.5, linewidth = 1.5) +
    guides(
      color = guide_legend(title = "Integer Copy Number"),
      alpha = guide_legend(title = "Weight")
    )
  if (missing(chr) || length(chr) > 1) {
    plot <- plot + facet_wrap(vars(chromosome),
      strip.position = "bottom", scales = "free_x",
      nrow = 1
    ) +
      theme_bw() +
      theme(
        axis.line = element_line(colour = "grey"),
        strip.placement = "outside",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white")
      )
  }
  plot
}

plot_one <- plot_cnvkit(cnr, cns, 19)
plot <- plot_cnvkit(cnr, cns, c("1", "2", "3"))
plot
