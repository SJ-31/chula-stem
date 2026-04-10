library(tidygraph)
library(igraph)
library(tidyverse)
library(checkmate)
library(ggraph)
library(pointblank)
library(here)

tb <- read_tsv(
  "/home/shannc/Bio_SDD/stem_synology/chula_mount/shannc/output/PHCase/PHcase_13/annotations/7-PHcase_13-VEP_small.tsv",
  col_types = list(Loc = "character")
) |>
  distinct(Loc, Alt, Ref, .keep_all = TRUE) |>
  separate_wider_delim(cols = "Loc", delim = ":", names = c("chr", "pos")) |>
  mutate(pos = as.numeric(pos))

gene_ref <- read_csv(here("analyses", "data", "ensembl_gene_data.csv"))

wanted <- "KRAS"

len <- 25250936 - 25205246

## cur <- tb |> filter(SYMBOL == )

#' Generate a chromosome graph to plot breakpoints
#'
#' @description
#' @param n_bins Number of nodes displayed for each chromosome. Incompatible with
#' `interval_length` parameter
make_gene_graph <- function(
  tb,
  gene_start,
  gene_stop,
  gene_col = "SYMBOL",
  hgsvg_col = "HGVSg",
  hgsvp_col = "HGVSp",
  hgsvc_col = "HGVSc",
  bin_width = 10^3,
  n_bins = NULL
) {
  assert_data_frame(tb, min.rows = 1)
  expect_col_vals_regex(tb, alt_col, regex = "[\\[\\]]")
  expect_col_vals_regex(tb, alt_col, regex = ":")

  tb <- rename(
    tb,
    gene = gene_col,
    hgsvg = hgsg_col,
    hgsvp = hgsp_col,
    hgsvc = hgsc_col
  ) |>
    distinct(hgsg, .keep_all = TRUE) |>
    mutate()

  # [2026-03-14 Sat] TODO: modify all of the below
  gene_length <- gene_stop - gene_start

  if (!is.null(n_bins)) {
    nodes <- tibble(
      pos = seq(0, gene_length, length.out = n_bins),
      id = seq(n_bins)
    )
  } else {
    nodes <- tibble(pos = seq(0, gene_length, by = bin_width)) |>
      mutate(id = seq_along(pos))
  }

  ## var_nodes <- tb |> select(gene, hgnc, hgsp, hgnc) |> mutate()

  name2key <- select(nodes, name, id) |> deframe()
  bp_edges <- tb |>
    mutate(
      from = paste0(chr, ":", pos),
      to = str_extract(alt, "([1-9]+:[0-9]+)", group = 1),
      ## direction = case_when(str_match(alt, "[") ~ "left" ),
    ) |>
    mutate(
      from = name2key[from],
      to = name2key[to],
      type = "breakpoint",
      label = alt
    )

  edges <- nodes |>
    group_by(chr) |>
    summarize(from = list(id), to = list(id)) |>
    mutate(
      from = lapply(from, \(x) head(x, -1)),
      to = lapply(to, \(x) tail(x, -1))
    ) |>
    unnest(c(from, to)) |>
    mutate(type = "scaffold") |>
    bind_rows(bp_edges) |>
    mutate(
      label = replace_values(label, NA ~ ""),
      color = case_when(type == "scaffold" ~ chr, .default = NA)
    )

  G <- tbl_graph(nodes = nodes, edges = edges)
}

## G <- make_chr_graph(breakpoints, bin_width = 10^6.8)

## chr2x <- local({
##   chr <- unique(vertex_attr(G)$chr)
##   setNames(seq(1, 24, length.out = length(chr)), chr)
## })

## x_y <- activate(G, nodes) |>
##   as_tibble() |>
##   mutate(x = chr2x[chr]) |>
##   group_by(chr) |>
##   group_split() |>
##   lapply(\(x) mutate(x, y = pos)) |>
##   bind_rows() |>
##   select(x, y)

## plot <- ggraph(G, "manual", x = x_y$x, y = x_y$y) +
##   geom_edge_link(
##     aes(color = color, label = label),
##     angle_calc = "along",
##     edge_width = 1.5,
##     label_dodge = unit(4, "mm"),
##     check_overlap = TRUE,
##     # TODO: add the correct arrows indicating directionality
##   ) +
##   geom_node_text(
##     aes(label = label),
##     nudge_x = -0.1,
##     hjust = "right",
##     fontface = ifelse(vertex_attr(G, "type") == "chr_label", "bold", "plain")
##   ) +
##   geom_node_point(size = vertex_attr(G)$size) +
##   guides(color = "none") + # TODO: this isn't working for some reason
##   theme_void() +
##   scale_y_reverse()
## plot
## ggsave(here(".scratch", "foobar.pdf"), plot, width = 10)
