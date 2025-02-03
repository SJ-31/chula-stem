library(tidyverse)
library(glue)
library(ggpubr)
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

plot_cnvkit_gene <- function(cnv_tb, locus) {

}

parse_chr_range <- function(chr_range) {
  splits <- str_split_1(chr_range, ":")
  chr <- splits[1]
  if (length(splits) > 1) {
    endpoints <- str_split_1(splits[2], "-")
    list(chr = chr, ir = IRanges(as.numeric(endpoints[1]), as.numeric(endpoints[2])))
  } else {
    list(chr = chr, ir = NULL)
  }
}

plot_cnvkit <- function(cns, cnr, chr, sizing, output) {
  cn_calls <- read_tsv(cns) |> format_chr_data()
  ratios <- read_tsv(cnr) |> format_chr_data()
  range_given <- FALSE
  if (!is.null(chr) && length(chr) >= 1) {
    if (str_detect(chr, "-")) range_given <- TRUE
    range <- parse_chr_range(chr)
    cn_calls <- cn_calls |> dplyr::filter(chromosome %in% range$chr)
    ratios <- ratios |> dplyr::filter(chromosome %in% range$chr)
  }
  merged <- ratios |> inner_join(cn_calls,
    by = join_by(chromosome, between(mid, y$start, y$end)), suffix = c("", "_right")
  )
  if (range_given) {
    vals <- map2(merged$start_right, merged$end_right, \(s, e) {
      intersect <- pintersect(IRanges(s, e), range$ir, resolve.empty = "start.x")
      if (width(intersect) > 0) {
        tibble(start_right = start(intersect), end_right = end(intersect))
      } else {
        tibble(start_right = NA, end_right = NA)
      }
    }) |> bind_rows()
    merged <- merged |>
      dplyr::select(-start_right, -end_right) |>
      bind_cols(vals) |>
      dplyr::filter(start >= start(range$ir) & end <= end(range$ir))
  }
  plot <- merged |>
    ggplot(aes(x = mid, y = log2, alpha = weight), stat = "identity") +
    ylab("Copy Number Ratio (log2)") +
    geom_point(size = 0.5) +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 15),
      minor_breaks = scales::minor_breaks_n(5),
      limits = c(min(ratios$log2), max(ratios$log2)),
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
      ncol = lget(sizing, "ncol", 3)
    ) +
      theme_bw() +
      theme(
        text = element_text(size = 15),
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
    without_chr <- str_replace(chr, "^chr", "")
    plot <- plot + xlab("Position (mb)") + labs(title = glue("Chromosome {without_chr}"))
  }
  if (range_given) {
    plot <- plot + scale_x_continuous(limits = c(range$start, range$end))
  }
  if (!is.null(output)) {
    ggsave(output, plot,
      dpi = lget(sizing, "dpi", 300), width = lget(sizing, "width", NA),
      height = lget(sizing, "height", NA), units = lget(sizing, "units", "in")
    )
  } else {
    plot
  }
}

pca_dgelist <- function(dgelist, plot_aes = list(shape = "group", color = "type"), log = TRUE) {
  if (class(dgelist) == "DGEList") {
    counts <- edgeR::cpm(dgelist$counts, log = log) # Convert counts into counts per million
    # When dgelist is normalized with `normLibSizes`, then cpm automatically generates
    # normalized counts from the computed normalization factors
    # This is roughly equivalent to the following:
    # 1. Get effective library size
    ## effective_lib_size <- dgelist$samples$lib.size * dgelist$samples$norm.factors
    # 2. Scale the counts by the lib size
    ## counts <- t(t(counts) / effective_lib_size)
  } else {
    counts <- dgelist$counts
  }
  pca <- prcomp(t(counts))
  summary <- summary(pca)
  var_pca1 <- summary$importance[, 1]["Proportion of Variance"]
  var_pca2 <- summary$importance[, 2]["Proportion of Variance"]
  pcs <- as_tibble(pca$x) |> bind_cols(dgelist$samples)
  plot_aes <- rlang::syms(plot_aes)
  aes_fn <- do.call(\(...) aes(x = PC1, y = PC2, ...), plot_aes)
  plot <- ggplot(pcs, aes_fn) +
    geom_point() +
    xlab(glue("PC1 ({var_pca1})")) +
    ylab(glue("PC2 ({var_pca2})"))
  plot
}
