library(here)
library(tidyverse)
library(glue)
## library(ensembldb)

cur_dir <- here("analyses", "brca")
env <- yaml::read_yaml(here(cur_dir, "env.yaml"))
mpath <- env$metadata_path

lget <- function(lst, key, default = NULL) {
  if (is.null(lst[[key]])) {
    default
  } else {
    lst[[key]]
  }
}

# Sample columns to keep or generate
shared_cols <- c(
  "join_d", # identifier column with which to join columns in raw data
  "subtype", # brca subtype
  "age",
  "treatment",
  "treatment_response",
  "platform",
  "histological_grade",
  "histological_type",
  "t_stage",
  "sample_type", # recurrent or primary. Default "primary"
  # Columns for markers
  "er_positive",
  "pr_positive"
)
# If the above columns cannot be found in `meta_remap` or `meta_fill`, they will be filled with NA

# Aggregate metadata
mdata <- lapply(names(env$datasets), \(name) {
  meta_files <- lget(env$datasets[[name]], "files", list())$meta
  if (!is.null(meta_files)) {
    tb <- lapply(meta_files, \(f) read_tsv(here(mpath, f))) |> bind_rows()
  } else {
    tb <- read_tsv(here(mpath, glue("{name}.tsv")))
  }
  join_id <- env$datasets[[name]]$id_col
  if (is.null(join_id)) {
    join_id <- env$id_col
  }
  tb <- tb |> mutate(dataset = name, join_id = !!as.symbol(join_id))
  if ("title" %in% colnames(tb) && class(tb$title) == "numeric") {
    tb$title <- as.character(tb$title)
  }
  tb
}) |>
  bind_rows()
