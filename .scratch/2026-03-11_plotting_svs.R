library(tidygraph)
library(igraph)
library(tidyverse)
library(checkmate)
library(ggraph)
library(pointblank)
library(here)

tb <- read_tsv(
  "/home/shannc/Bio_SDD/stem_synology/chula_mount/shannc/output/PHCase/PHcase_21/annotations/8-PHcase_21-VEP_SV.tsv",
  col_types = list(Loc = "character")
) |>
  distinct(Loc, Alt, Ref, .keep_all = TRUE) |>
  separate_wider_delim(cols = "Loc", delim = ":", names = c("chr", "pos")) |>
  mutate(pos = as.numeric(pos))

CHR_LENGTHS <- c(
  "1" = 248956422,
  "2" = 242193529,
  "3" = 198295559,
  "4" = 190214555,
  "5" = 181538259,
  "6" = 170805979,
  "7" = 159345973,
  "8" = 145138636,
  "9" = 138394717,
  "10" = 133797422,
  "11" = 135086622,
  "12" = 133275309,
  "13" = 114364328,
  "14" = 107043718,
  "15" = 101991189,
  "16" = 90338345,
  "17" = 83257441,
  "18" = 80373285,
  "19" = 58617616,
  "20" = 64444167,
  "21" = 46709983,
  "22" = 50818468,
  "X" = 156040895,
  "Y" = 57227415
)

breakpoints <- tb |>
  filter(SVTYPE == "BND") |>
  select(chr, pos, Ref, Alt)


#' Generate a chromosome graph to plot breakpoints
#'
#' @description
#' @param n_bins Number of nodes displayed for each chromosome. Incompatible with
#' `interval_length` parameter
make_chr_graph <- function(
  tb,
  chr_col = "chr",
  alt_col = "Alt",
  bin_width = 10^6,
  n_bins = NULL
) {
  assert_data_frame(tb, min.rows = 1)
  expect_col_vals_regex(tb, alt_col, regex = "[\\[\\]]")
  expect_col_vals_regex(tb, alt_col, regex = ":")

  tb <- rename(tb, chr = chr_col, alt = alt_col) |>
    distinct(chr, pos, .keep_all = TRUE)

  chrs_present <- c(
    unique(tb$chr),
    str_extract(tb$alt, "([1-9]+):", group = 1)
  )
  kept_chr <- CHR_LENGTHS[names(CHR_LENGTHS) %in% chrs_present]

  bp_nodes <- bind_rows(
    select(tb, chr, pos),
    select(
      mutate(
        select(tb, alt),
        chr = str_extract(alt, "([1-9]+):", group = 1),
        pos = as.numeric(str_extract(alt, ":([0-9]+)", group = 1))
      ),
      -alt
    )
  ) |>
    mutate(type = "breakpoint", name = paste0(chr, ":", pos))

  if (!is.null(n_bins)) {
    nested <- mutate(
      enframe(kept_chr, name = "chr", value = "length"),
      pos = lapply(length, \(l) seq(0, l, length.out = n_bins)),
      seq_id = lapply(length, \(.) seq(n_bins))
    )
  } else {
    nested <- mutate(
      enframe(kept_chr, name = "chr", value = "length"),
      pos = lapply(length, \(l) seq(0, l, by = bin_width)),
      seq_id = lapply(pos, \(p) seq_along(p))
    )
  }

  nodes <- nested |>
    unnest(c(pos, seq_id)) |>
    mutate(name = paste0(chr, ".", seq_id), type = "scaffold") |>
    select(-seq_id) |>
    bind_rows(bp_nodes) |>
    group_by(chr) |>
    group_split() |>
    lapply(\(x) {
      cur_chr <- head(x$chr, n = 1)
      add_row(arrange(x, pos), chr = cur_chr, pos = -1000, type = "chr_label")
    }) |>
    bind_rows() |>
    mutate(id = row_number()) |>
    mutate(
      label = case_when(
        pos == 0 ~ "",
        type == "scaffold" ~ as.character(round(pos / 10^6, 0)),
        type == "chr_label" ~ paste0("chr", chr),
        .default = ""
      ),
      size = case_when(
        type == "breakpoint" ~ 1,
        type == "chr_label" ~ 2,
        .default = 0,
      )
    )

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


# [2026-03-13 Fri] TODO: Getting there, but you should experiment with other layouts too
# Also get rid of the nodes that aren't breakpoints

G <- make_chr_graph(breakpoints, bin_width = 10^6.8)

chr2x <- local({
  chr <- unique(vertex_attr(G)$chr)
  setNames(seq(1, 24, length.out = length(chr)), chr)
})

x_y <- activate(G, nodes) |>
  as_tibble() |>
  mutate(x = chr2x[chr]) |>
  group_by(chr) |>
  group_split() |>
  lapply(\(x) mutate(x, y = pos)) |>
  bind_rows() |>
  select(x, y)

plot <- ggraph(G, "manual", x = x_y$x, y = x_y$y) +
  geom_edge_link(
    aes(color = color, label = label),
    angle_calc = "along",
    edge_width = 1.5,
    label_dodge = unit(4, "mm"),
    check_overlap = TRUE,
    # TODO: add the correct arrows indicating directionality
  ) +
  geom_node_text(
    aes(label = label),
    nudge_x = -0.1,
    hjust = "right",
    fontface = ifelse(vertex_attr(G, "type") == "chr_label", "bold", "plain")
  ) +
  geom_node_point(size = vertex_attr(G)$size) +
  guides(color = "none") + # TODO: this isn't working for some reason
  theme_void() +
  scale_y_reverse()

plot
ggsave(here(".scratch", "foobar.pdf"), plot, width = 10)


chr_line_plot <- function(
  tb,
  chr_col = "chr",
  alt_col = "Alt",
  bin_width = 10^6,
  n_bins = NULL
) {
  assert_data_frame(tb, min.rows = 1)
  expect_col_vals_regex(tb, alt_col, regex = "[\\[\\]]")
  expect_col_vals_regex(tb, alt_col, regex = ":")

  tb <- rename(tb, chr = chr_col, alt = alt_col) |>
    distinct(chr, pos, .keep_all = TRUE)

  chrs_present <- c(
    unique(tb$chr),
    str_extract(tb$alt, "([1-9]+):", group = 1)
  )
  kept_chr <- CHR_LENGTHS[names(CHR_LENGTHS) %in% chrs_present]

  chr_x_map <- local({
    lvls <- levels(as.factor(names(kept_chr)))
    seq_along(lvls) |> `names<-`(lvls)
  })

  chr_table <- enframe(kept_chr, name = "chr", value = "y")
  chr_table <- bind_rows(chr_table, mutate(chr_table, y = 0))

  bp_table <- select(tb, chr, pos, alt) |>
    mutate(
      chr_to = str_extract(alt, "([0-9]+):", group = 1),
      pos = as.numeric(pos),
      chr_to_pos = as.numeric(str_extract(alt, ":([0-9]+)", group = 1)),
      from_x = chr_x_map[chr],
      to_x = chr_x_map[chr_to],
      alt = str_replace(alt, ":([0-9]+)", replacement = \(x) {
        val <- str_remove(x, ":") |> as.numeric()
        paste0(":", formatC(val, digits = 2, format = "e"))
      })
    )

  ggplot(chr_table, aes(x = chr, y = y, color = chr)) +
    geom_line(linewidth = 3) +
    geom_segment(
      data = bp_table,
      aes(x = from_x, xend = to_x, y = pos, yend = chr_to_pos)
    ) +
    geom_label(
      data = bp_table,
      aes(x = from_x, label = alt, y = pos),
      color = "black",
      fill = "white",
      hjust = "right"
    ) +
    theme_void() +
    theme(
      axis.text.y = element_text(),
      axis.title.y = element_text(face = "bold", angle = 90),
    ) +
    scale_y_reverse() +
    guides(color = guide_legend(title = "Chromosome")) +
    ylab("Position")
  bp_table
}
# %%

library(circlize)

sv_circos <- function(
  tb,
  chr_col = "chr",
  alt_col = "Alt",
  genome = "hg38"
) {
  tb <- rename(tb, chr = chr_col, alt = alt_col) |>
    distinct(chr, pos, .keep_all = TRUE)

  chrs_present <- c(
    unique(tb$chr),
    str_extract(tb$alt, "([0-9]+):", group = 1)
  ) |>
    unique()

  bp_table <- select(tb, chr, pos, alt) |>
    mutate(
      chr = if_else(
        str_starts(chr, "chr"),
        chr,
        paste0("chr", chr)
      ),
      pos = as.numeric(pos),
      chr_to = str_extract(alt, "([0-9]+):", group = 1),
      chr_to = if_else(
        str_starts(chr_to, "chr"),
        chr_to,
        paste0("chr", chr_to)
      ),
      pos = as.numeric(pos),
      chr_to_pos = as.numeric(str_extract(alt, ":([0-9]+)", group = 1)),
      alt = str_replace(alt, ":([0-9]+)", replacement = \(x) {
        val <- str_remove(x, ":") |> as.numeric()
        paste0(":", formatC(val, digits = 2, format = "e"))
      })
    )

  to <- bp_table |> select(chr, pos, pos)
  from <- bp_table |> select(chr_to, chr_to_pos, chr_to_pos)

  circos.initializeWithIdeogram(
    species = genome,
    chromosome.index = paste0("chr", chrs_present)
  )
  circos.genomicLink(to, from)
}

# [2026-03-26 Thu] it's alright, but better to use something that displays the joinings
# clearer
sv_circos(breakpoints)
