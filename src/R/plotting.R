library(tidyverse)
library(glue)
library(grid)
library(ggplot2)

#' Returns "default" if the value of "key" in "list" is NULL
#'
lget <- function(list, key, default) {
  val <- list[[key]]
  if (is.null(val)) {
    default
  } else {
    val
  }
}

format_chr_data <- function(tb, divide_by = 1000) {
  chr_levels <- c(as.character(1:22), "X", "Y") %>% paste0("chr", .)
  if (!str_detect(tb$chromosome[1], "^chr")) {
    tb$chromosome <- paste0("chr", tb$chromosome)
  }
  tb |> mutate(
    start = start / divide_by,
    end = end / divide_by,
    mid = (end - start) / 2 + start,
    chromosome = factor(tb$chromosome, levels = chr_levels)
  )
}

plot_cnvkit <- function(cns, cnr, chr, sizing, output) {
  cn_calls <- read_tsv(cns) |> format_chr_data()
  ratios <- read_tsv(cnr) |> format_chr_data()
  if (!is.null(chr) && length(chr) >= 1) {
    cn_calls <- cn_calls |> filter(chromosome %in% chr)
    ratios <- ratios |> filter(chromosome %in% chr)
  }
  plot <- ratios |>
    inner_join(cn_calls,
      by = join_by(chromosome, between(mid, y$start, y$end)),
      suffix = c("", "_right")
    ) |>
    ggplot(aes(x = mid, y = log2, alpha = weight), stat = "identity") +
    ylab("Copy Number Ratio (log2)") +
    geom_point() +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 15),
      minor_breaks = scales::minor_breaks_n(5)
    ) +
    geom_segment(aes(
      x = start_right, xend = end_right,
      y = log2_right, yend = log2_right,
      color = factor(cn)
    ), alpha = 0.5, linewidth = 1.5) +
    guides(
      color = guide_legend(title = "Integer Copy Number"),
      alpha = guide_legend(title = "Weight")
    ) +
    xlab(NULL)
  if (is.null(chr) || length(chr) > 1) {
    without_chr <- str_replace(chr, "^chr", "")
    title <- ifelse(missing(chr), "All Chromosomes",
      glue("Chromosomes {paste(without_chr, collapse = ', ')}")
    )
    plot <- plot + facet_wrap(vars(chromosome),
      strip.position = "bottom", scales = "free_x",
      nrow = 1
    ) +
      theme_bw() +
      theme(
        axis.line = element_line(colour = "grey"),
        strip.placement = "outside",
        panel.grid.minor.x = element_blank(),
        legend.title = element_text(face = "bold"),
        title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        axis.ticks.y = element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white")
      ) + labs(title = title)
  } else {
    plot <- plot + xlab("Position (mb)") + labs(title = glue("Chromosome {chr}"))
  }
  ggsave(output, plot,
    dpi = lget(sizing, "dpi", 300), width = lget(sizing, "width", NA),
    height = lget(sizing, "height", NA), units = lget(sizing, "units", "in")
  )
}
