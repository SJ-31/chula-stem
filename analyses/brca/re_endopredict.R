suppressMessages({
  library(here)
  library(tidyverse)
  library(glue)
  library(AnnotationHub)
  library(Biobase)
  library(GEOquery)
  setAnnotationHubOption("CACHE", here(".cache", "AnnotationHub"))
})
geo_cache <- here(".cache", "GEO")
date <- format(Sys.time(), "%Y-%m-%d")
cur_dir <- here("analyses", "brca")
env <- yaml::read_yaml(here(cur_dir, "env.yaml"))
outdir <- here("analyses", "output", "brca", "re_endopredict")
dir.create(outdir)
mpath <- env$metadata_path

# Authors already re-normalized the data from the different cohorts in the same way as theirs
# shared variable is time to distant recurrence

# Unit of dmfs_time is months
shared_cols <- c("geo_accession", "dmfs_time", "dmfs_status")

mfiles <- list(
  GSE26971 = "GSE26971.tsv",
  GSE12093 = "GSE12093.tsv",
  `GSE6532-GPL96` = "GSE6532_LUMINAL_demo.txt"
)

efiles <- c(
  "GSE26971_GSE12093_dataset.txt",
  "GSE26971_GSE6532_dataset.txt",
  "GSE26971.tsv"
)

expr_dir <- here(env$raw_path, "GSE26971")

## * Gather metadata
meta <- lapply(names(mfiles), \(x) {
  file <- mfiles[[x]]
  tb <- read_tsv(here(mpath, file))
  spec <- env$datasets[[x]]
  to_remap <- unlist(spec$meta_remap)
  if (!is.null(spec$id_col)) {
    tb <- dplyr::rename(tb, geo_accession = spec$id_col)
  }
  if (!is.null(to_remap)) tb <- dplyr::rename(tb, to_remap)
  if (x == "GSE6532-GPL96") {
    tb <- mutate(tb, dmfs_time = dmfs_time / (365.24 / 12))
  }
  to_replace <- spec$meta_replace
  for (col in names(to_replace)) {
    if (is.character(tb[[col]])) {
      vec <- str_to_lower(tb[[col]])
    } else {
      vec <- tb[[col]]
    }
    tb[[col]] <- unlist(to_replace[[col]])[vec] |>
      unlist(use.names = FALSE)
    if (!is.null(to_replace[[".default"]])) {
      tb[[col]] <- replace_na(tb[[col]], to_replace[[".default"]])
    }
  }
  if (!"subcohort" %in% to_remap) tb <- mutate(tb, subcohort = x)
  tb |>
    mutate(tb, dmfs_time = floor(dmfs_time)) |>
    dplyr::select(all_of(shared_cols))
}) |>
  bind_rows()

## * Gather expression

expr_tmp <- lapply(efiles, \(x) read_tsv(here(expr_dir, x))) |>
  reduce(\(x, y) left_join(x, y, by = join_by("ID_REF")))
