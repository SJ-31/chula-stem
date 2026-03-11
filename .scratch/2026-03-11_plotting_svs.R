library(tidygraph)
library(tidyverse)
library(ggraph)
library(here)

tb <- read_tsv(
  "/home/shannc/Bio_SDD/stem_synology/chula_mount/shannc/output/PHCase/PHcase_13/annotations/8-PHcase_13-VEP_SV.tsv",
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

n_bins <- 15

nodes <- enframe(CHR_LENGTHS, name = "chr", value = "length") |>
  mutate(
    segments = lapply(length, \(l) seq(0, l, length.out = n_bins)),
    seq_id = lapply(length, \(.) seq(n_bins))
  ) |>
  unnest(c(segments, seq_id)) |>
  mutate(name = paste0(chr, ".", seq_id)) |>
  select(-seq_id)

# TODO: this works well, now just add in the breakpoints before making the edges

edges <- nodes |>
  mutate(id = row_number()) |>
  group_by(chr) |>
  summarize(from = list(id), to = list(id)) |>
  mutate(
    from = lapply(from, \(x) head(x, -1)),
    to = lapply(to, \(x) tail(x, -1))
  ) |>
  unnest(c(from, to))

G <- tbl_graph(nodes = nodes, edges = edges)


ggraph(G) + geom_node_point() + geom_edge_fan(aes(color = chr)) + theme_void()

## test_g <- tbl_graph(
##   nodes = tibble(
##     name = c(
##       "11:72109776",
##       "3:105649604",
##       "6:156547959",
##       "11:START",
##       "11:END",
##       "3:START",
##       "3:END"
##       "6:START"
##     )
##   ),
##   edges = tibble(from = )
## )
