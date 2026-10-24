observed_de <- read_csv(
  glue("{env$data_root}/cohort/annotations/samples_de.csv")
) |>
  filter(proba_de > 0.6) |>
  rename(lfc = lfc_median)

de_tb <- observed_de |>
  filter(analysis_group == "scVI_subtype") |>
  distinct(gene, .keep_all = TRUE)
# Just for demonstratio

cur_g <- activate(fis_g, nodes) |>
  left_join(de_tb, by = join_by(x$name == y$gene)) |>
  mutate(is_de = !is.na(proba_de)) |>
  keep_interesting_comps("is_de", 2) |>
  purrr::reduce(\(x, y) bind_graphs(x, y))

plot_de_graph(cur_g)
