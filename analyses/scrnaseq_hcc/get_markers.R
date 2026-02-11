library(tidyverse)
library(glue)
library(here)

cell_markers <- list()

outdir <- here("analyses", "data", "hcc_reference")

## * Markers from previous project
wanted_cellmarker2 <- list(
  "t_helper_17_cell",
  "gamma_delta_t_cell",
  "m1_macrophage",
  "t_helper_1_cell",
  "m2_macrophage",
  "exhausted_t_cell",
  "stromal_cell",
  "natural_killer_cell",
  "cancer_cell",
  "neutrophil",
  "exhausted_cd8+_t_cell",
  "macrophage",
  "monocyte",
  "hepatocyte"
)
# Set up path
too_predict <- file.path(dirname(here()), "too-predict")

# Read and prepare ID mapping
id_map <- read_tsv(
  file.path(
    too_predict,
    "data/reference/Homo_sapiens.GRCh38.113.gene_id_mapping.tsv"
  )
) %>%
  rename(hgnc_symbol = gene_name)

# Read cell markers and gene sets
cell_markers_old <- read_tsv(
  file.path(too_predict, "data/reference/cell_markers_custom_meta.tsv")
)

gene_sets_old <- yaml::read_yaml(
  file.path(too_predict, "data/reference/gene_sets_custom.yaml")
)

# Initialize cell_markers list

joined_data <- cell_markers_old %>%
  left_join(id_map, by = c("ensembl" = "gene_id")) %>%
  filter(!is.na(hgnc_symbol))

for (cell_type in unique(joined_data$cell_type)) {
  if (cell_type %in% wanted_cellmarker2) {
    genes <- joined_data %>%
      filter(cell_type == !!cell_type) %>%
      pull(hgnc_symbol) %>%
      unique()
    cell_markers[[glue("cellmarker2-{cell_type}")]] <- genes
  }
}

gene_sets <- enframe(gene_sets_old, name = "set_name", value = "gene_id") |>
  unnest(gene_id) |>
  left_join(id_map, by = join_by(gene_id)) |>
  filter(!is.na(hgnc_symbol)) |>
  group_by(set_name) |>
  summarise(hgnc_symbol = list(hgnc_symbol)) |>
  filter(length(hgnc_symbol) >= 3) |>
  deframe()

## * DISCO

# By definition, the full list of markers for a child cell type includes the markers of its parents
disco <- read_csv(here("analyses", "data", "disco_cell_types.csv")) |>
  rename_with(
    \(x) str_remove_all(str_to_lower(str_replace_all(x, " ", "_")), "[()]")
  ) |>
  mutate(
    marker_curated = lapply(
      marker_curated,
      \(x) if (!is.na(x)) str_split_1(x, ";") else NA
    )
  )

get_parent_markers <- function(child_type, disco_tb) {
  rec_helper <- function(ct, markers) {
    tb <- disco_tb |> filter(cell_type == ct)
    markers <- c(markers, tb$marker_curated)
    if (tb$parent == "Cell Type" || nrow(tb) == 0) {
      markers
    } else {
      rec_helper(tb$parent, markers)
    }
  }
  rec_helper(child_type, c()) |>
    unlist() |>
    discard(is.na)
}

disco_hepatocytes <- disco |>
  filter(grepl("[Hh]epatocyte", cell_type) & cell_type != "Hepatocyte")
disco_cholangiocytes <- disco |> filter(grepl("[cC]holangiocyte", cell_type))

for (ct in disco_cholangiocytes$cell_type) {
  cell_markers[[glue("disco-{ct}")]] <- unique(get_parent_markers(ct, disco))
}

for (ct in disco_hepatocytes$cell_type) {
  cell_markers[[glue("disco-{ct}")]] <- get_parent_markers(ct, disco)
}

## * MSIGDB
# %%

msigdb <- msigdbr::msigdbr() |> filter(gs_source_species == "HS")

from_liver <- msigdb |>
  filter(grepl("AIZARANI_LIVER.*", gs_name)) |>
  distinct(gene_symbol, .keep_all = TRUE)

liver_cell2rx <- list(
  hepatocyte = "HEPATOCYTES",
  kupffer_cell = "KUPFFER_CELLS",
  stellate_cell = "STELLATE_CELLS",
  epcam_pos_cholangiocyte = "EPCAM_POS_BILE_DUCT",
  sinusoidal_endothelial_cell = "LSECS",
  liver_rsident_nkt = "NK_NKT_CELLS",
  liver_resident_b_cell = "RESIDENT_B_CELLS"
)
for (n in names(liver_cell2rx)) {
  cell_markers[[glue("msigdb_aizarani-{n}")]] <- from_liver |>
    filter(grepl(liver_cell2rx[[n]], gs_name)) |>
    pluck("gene_symbol") |>
    unique()
}

## * Write reference files

yaml::write_yaml(cell_markers, here(outdir, "cell_markers.yaml"))
yaml::write_yaml(gene_sets, here(outdir, "gene_sets.yaml"))
