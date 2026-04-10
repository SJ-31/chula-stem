library(tidyverse)
library(igraph)
library(visNetwork)
library(reticulate)
library(checkmate)
library(here)
library(tidygraph)
use_condaenv("stem-base")
source(here("src", "R", "utils.R"))

source(here("src", "R", "plotting.R"))

workdir <- here("analyses", "scrnaseq_hcc")
fn <- new.env()

source_python(here(workdir, "functions.py"), envir = fn)
ad <- import("anndata")
sc <- import("scanpy")
dc <- import("decoupler")
sc_utils <- import("chula_stem.sc_rnaseq")
data <- here("analyses", "data")
# %%

sys <- import("sys")
sys$argv <- c(sys$argv, "test=True")
yte <- import("yte")
env <- yte$process_yaml(paste0(
  read_lines(here(workdir, "cellranger_config.yaml")),
  collapse = "\n"
))


go_term2id <- read_csv(env$go_term2id) |> distinct(term, .keep_all = TRUE)
go <- read_graph(env$go_graph, "gml") |>
  as_tbl_graph() |>
  select(-id) |>
  rename(distance_to_ns = "distancetons")

## ** Reactome
observed <- read_csv(
  glue("{env$data_root}/cohort/annotations/clusters_gene_set_activity.csv")
) |>
  filter(str_starts(name, "Reactome:")) |>
  mutate(name = str_remove(name, "^Reactome:"))

observed_go <- read_csv(
  glue("{env$data_root}/cohort/annotations/clusters_gene_set_activity.csv")
) |>
  filter(str_starts(name, "GO_..:")) |>
  mutate(name = str_remove(name, "GO_..:")) |>
  inner_join(go_term2id, by = join_by(x$name == y$term))


cur <- "leiden_res0.1"
cg <- 4

tb <- observed |> filter(clustering == cur & group == cg)

tb2 <- observed_go |>
  filter(clustering == cur & group == cg) |>
  select(-name) |>
  rename(name = "id") |>
  distinct(name, .keep_all = TRUE)

# TODO: loop over each of the components below, write a pdf file,
# and merge the pdfs together for each grouping of clustering and group

# TODO: can you rename the pages?

rg <- prepare_reactome_pwy(
  "/home/shannc/Bio_SDD/chula-stem/analyses/data/ReactomePathways.txt",
  "/home/shannc/Bio_SDD/chula-stem/analyses/data/ReactomePathwaysRelation.txt"
) |>
  activate(nodes) |>
  left_join(tb, by = join_by(name)) |>
  mutate(enriched = !is.na(padj))

rg_e <- keep_interesting_comps(rg, "enriched")


# Prune nodes until all leaf nodes are enriched
max_length <- 20
tmp <- rg_e[[6]] |>
  keep_interesting_leaves("enriched") |>
  activate(nodes) |>
  mutate(
    name = str_wrap(name, width = max_length)
  )


plot_enriched_graph(tmp)

ggsave(here(workdir, "graph.pdf"), width = 12, height = 12)

## ** De genes again

obs_de <- read_csv(glue(
  "{env$data_root}/cohort/annotations/clusters-scVI_de.csv"
))

fi_g <- prepare_reactome_fi(env$reactome_fi)

cur_clust <- obs_de |>
  filter(
    contrast == "1 vs Rest" &
      clustering == "leiden_res0.1" &
      proba_de > 0.9 &
      abs(lfc_mean) > 1.5
  ) |>
  rename(lfc = "lfc_median")


## * Interactive visualization
# %%

vnodes <- as_tibble(tmp, "nodes") |>
  mutate(
    id = row_number(),
    color = recode_values(enriched, FALSE ~ "grey", TRUE ~ "red")
  ) |>
  rename(tooltip = "name") |>
  as.data.frame()
edges <- as_tibble(tmp, "edges") |> as.data.frame()

visNetwork(vnodes, edges) |>
  visPhysics(stabilization = FALSE) |>
  visEdges(smooth = FALSE)
