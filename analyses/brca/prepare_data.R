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

##' Unify representation of marker gene status into binary vector
recode_status <- function(vec) {
  case_match(
    str_to_lower(as.character(vec)),
    "positive" ~ 1,
    "negative" ~ 0,
    "neg" ~ 0,
    "pos" ~ 1,
    "yes" ~ 1,
    "no" ~ 0,
    "p" ~ 1,
    "n" ~ 0,
    .default = 0
  )
}

##' Unify subtype names
recode_subtype <- function(vec) {
  vec <- str_to_lower(vec) |> str_replace_all(" ", "_")
  case_match(
    vec,
    "luma" ~ "luminal_a",
    "lumb" ~ "luminal_b",
    "lum" ~ "luminal",
    "her_2" ~ "her2",
    "tn" ~ "triple_negative",
    .default = vec
  )
}

# Sample columns to keep or generate
SHARED_COLS <- c(
  "join_id", # identifier column with which to join columns in raw data
  "dataset",
  "subtype", # brca subtype
  "age",
  "treatment", # "none" if no treament given
  "treatment_response",
  "platform",
  "histological_grade",
  "histological_type",
  "t_stage",
  "sample_type", # recurrent or primary. Default "primary"
  # Columns for markers
  "er_status",
  "pr_status",
  "her2_status"
)
# If the above columns cannot be found in `meta_remap` or `meta_fill`, they will be filled with NA

MARKER_COLS <- keep(SHARED_COLS, \(x) str_detect(x, "_status$"))
TO_CHARACTER <- c("histological_grade")

# Aggregate metadata
mdata <- lapply(names(env$datasets), \(name) {
  cur <- env$datasets[[name]]
  meta_files <- lget(cur, "files", list())$meta
  if (!is.null(meta_files)) {
    tb <- lapply(meta_files, \(f) suppressMessages(read_tsv(here(mpath, f)))) |>
      bind_rows()
  } else {
    tb <- suppressMessages(read_tsv(here(mpath, glue("{name}.tsv"))))
  }
  join_id <- cur$id_col
  if (is.null(join_id)) {
    join_id <- env$id_col
  }
  tb <- tb |> mutate(dataset = name, join_id = !!as.symbol(join_id))
  to_remap <- unlist(cur$meta_remap)
  to_fill <- cur$meta_fill
  if (!is.null(to_remap)) {
    tb <- dplyr::rename(tb, to_remap)
  }
  if (!is.null(to_fill)) {
    tb <- mutate(tb, !!!to_fill)
  }
  remaining_cols <- setdiff(SHARED_COLS, c(names(to_remap), names(to_fill)))
  add_na <- rep(NA, length(remaining_cols)) |>
    as.list() |>
    setNames(remaining_cols)

  to_replace <- cur$meta_replace
  for (col in names(to_replace)) {
    if (is.character(tb[[col]])) {
      vec <- str_to_lower(tb[[col]])
    } else {
      vec <- tb[[col]]
    }
    tb[[col]] <- unlist(to_replace[[col]])[vec] |>
      unlist(use.names = FALSE)
  }
  tb <- mutate(tb, !!!add_na)
  tb |>
    select(all_of(SHARED_COLS)) |>
    mutate(
      t_stage = as.character(t_stage),
      sample_type = as.character(sample_type),
      subtype = recode_subtype(subtype)
    ) |>
    mutate(across(all_of(MARKER_COLS), recode_status)) |>
    mutate(across(all_of(TO_CHARACTER), as.character))
}) |>
  bind_rows()
