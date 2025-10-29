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
outdir <- here("analyses", "output", "brca")
mpath <- env$metadata_path
# %%

## * Helper functions

# Extract expression data in sample tables of GSE80999
format_GSE80999 <- function() {
  file <- here(env$raw_path, "GSE80999", "GSE80999_family.soft.gz")
  gsms <- GSMList(getGEO(filename = file))
  lapply(names(gsms), \(name) {
    Table(gsms[[name]]) |>
      as_tibble() |>
      rename(VALUE = name)
  }) |>
    purrr::reduce(\(x, y) left_join(x, y, by = join_by(ID_REF))) |>
    write_tsv(here(env$raw_path, "GSE80999", "GSE80999.tsv"))
}

lget <- function(lst, key, default = NULL) {
  if (is.null(lst[[key]])) {
    default
  } else {
    lst[[key]]
  }
}


#' Replace all values in character vector `x` if they match a regexp in `spec`
#'
#' @param x character vector
#' @param spec list whose names are values to replace with and values are the corresponding regexps
replace_re_matches <- function(x, spec, .default = NA) {
  match_vecs <- lapply(names(spec), \(s) {
    matches <- str_detect(x, spec[[s]])
    vals <- rep(NA, length(matches))
    vals[matches] <- s
    vals
  })
  replace_na(coalesce(!!!match_vecs), .default)
}

recode_treatment_response <- function(vec) {
  vec <- str_to_lower(vec)
  case_match(
    vec,
    "pcr" ~ "pCR",
    "cr+kl" ~ "CR+known_lesion", # WARNING: this is an AI-assisted guess
    "cr" ~ "CR",
    "nocr" ~ "no_CR",
    "no pcr" ~ "no_pCR",
    "rd" ~ "no_CR",
    "pd" ~ "PD",
    "ne" ~ NA,
    .default = vec
  )
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
    .default = NA
  )
}

platform_ids2name <- function(ids) {
  ids <- as.character(ids)
  case_match(
    ids,
    "GPL20078" ~ "Agendia32627_DPv1.14_SCFGplus",
    "GPL6884" ~ "Illumina HumanWG-6 v3.0 expression beadchip",
    "GPL30493" ~ "Agilent_human_DiscoverPrint_15746",
    .default = ids
  )
}

recode_t_stage <- function(x) {
  # Keep prefixes
  x <- str_trim(x)
  case_when(
    str_detect(x, "^[1-4]$") ~ paste0("T", x),
    str_detect(x, "[()]") ~ str_remove_all(x, "\\(.*\\)"),
    x == "999" ~ NA, # Weird formatting in GSE123845
    .default = x
  )
}

recode_hist_type <- function(x) {
  x <- str_to_lower(x) |>
    str_remove_all(",") |>
    str_replace_all(" ", "_")
  case_match(x, NA ~ NA, .default = x) # TODO: might want to ask for help on this
}

recode_hist_grade <- function(x) {
  x <- str_trim(str_remove_all(x, "SBR"))
  case_match(x, "I" ~ "1", "II" ~ "2", "III" ~ "3", .default = x)
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

## * Metadata columns

# Sample columns to keep or generate
SHARED_COLS <- c(
  "join_id", # identifier column with which to join columns in raw data
  "patient_id", # identifier column for patients (some datasets have pre-, post-)
  "dataset",
  "subtype", # brca subtype
  "age",
  "treatment", # "none" if no treament given
  "treatment_response",
  "platform",
  "platform_name",
  "histological_grade",
  "histological_type",
  "collection_period",
  "t_stage",
  "sample_type",
  "pre_post",
  "recurrent", # binary marker for recurrent or not
  # Columns for molecular markers
  "er_status",
  "pr_status",
  "her2_status"
)
# If the above columns cannot be found in `meta_remap` or `meta_fill`, they will be filled with NA

MARKER_COLS <- keep(SHARED_COLS, \(x) str_detect(x, "_status$"))
TO_CHARACTER <- c("join_id", "patient_id", "collection_period")
TO_FACTOR <- c("histological_grade", "t_stage", "sample_type")

## * Aggregate metadata
mdata <- lapply(names(env$datasets), \(name) {
  cur <- env$datasets[[name]]
  meta_files <- lget(cur, "files", list())$meta
  if (!is.null(meta_files)) {
    tb <- lapply(meta_files, \(f) suppressMessages(read_tsv(here(mpath, f)))) |>
      bind_rows()
  } else {
    tb <- suppressMessages(read_tsv(here(mpath, glue("{name}.tsv"))))
  }
  to_remap <- unlist(cur$meta_remap)
  to_remap["join_id"] <- ifelse(!is.null(cur$id_col), cur$id_col, env$id_col)
  if (is.na(to_remap["platform"]) && env$platform_col %in% colnames(tb)) {
    to_remap["platform"] <- env$platform_col
  }
  if (!is.null(to_remap)) {
    tb <- dplyr::rename(tb, to_remap)
  }
  if (!"patient_id" %in% colnames(tb)) {
    tb <- mutate(tb, patient_id = join_id)
  }
  to_fill <- cur$meta_fill
  if (!is.null(to_fill)) {
    tb <- mutate(tb, !!!to_fill)
  }
  to_filter_out <- cur$meta_filter
  for (col in names(to_filter_out)) {
    blacklist <- to_filter_out[[col]]
    tb <- tb[!tb[[col]] %in% blacklist, ]
  }
  remaining_cols <- setdiff(
    SHARED_COLS,
    c(names(to_remap), names(to_fill), "join_id")
  )
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
    if (!is.null(to_replace[[".default"]])) {
      tb[[col]] <- replace_na(tb[[col]], to_replace[[".default"]])
    }
  }

  match_re <- cur$meta_match_re
  for (col in names(match_re)) {
    tb[[col]] <- replace_re_matches(tb[[col]], spec = match_re[[col]])
  }

  to_remove <- cur$meta_str_remove
  for (col in names(to_remove)) {
    tb[[col]] <- purrr::reduce(
      to_remove[[col]],
      \(acc, pat) str_remove_all(acc, pat),
      .init = tb[[col]]
    )
  }

  tb <- mutate(tb, !!!add_na) |>
    dplyr::select(all_of(SHARED_COLS)) |>
    mutate(
      subtype = recode_subtype(subtype),
      dataset = name,
      platform_name = platform_ids2name(platform),
      platform_type = if (!is.null(cur$microarray)) "microarray" else "RNA-seq",
      t_stage = recode_t_stage(t_stage),
      treatment_response = recode_treatment_response(treatment_response),
      histological_type = recode_hist_type(histological_type),
      histological_grade = recode_hist_grade(histological_grade),
      sample_type = replace_na(str_to_lower(sample_type), "primary"),
      recurrent = replace_na(recurrent, 0)
    ) |>
    mutate(across(all_of(MARKER_COLS), recode_status)) |>
    mutate(across(all_of(TO_CHARACTER), as.character)) |>
    mutate(across(all_of(TO_FACTOR), as.factor)) |>
    mutate(across(
      where(is.character) | where(is.factor),
      \(x) iconv(x, from = "UTF-8", to = "ASCII", sub = "")
    ))
  tb
}) |>
  bind_rows()
# %%

# TODO: need gene id mappings for
# and microarrays
# - "Agilent-028004 SurePrint G3 Human GE 8x60K Microarray (Probe Name Version)"
# - "Illumina HumanWG-6 v3.0 expression beadchip"
# - "Agilent_human_DiscoverPrint_15746"
# - Agendia32627_DPv1.14_SCFGplus

## * Aggregate samples

ah <- AnnotationHub(localHub = TRUE)
## ah <- AnnotationHub()
db <- query(ah, "org.Hs.eg.db")[[1]]
# %%

esets <- lapply(names(env$datasets), \(name) {
  cur <- env$datasets[[name]]
  gene_id_style <- cur$gene_id
  if (is.null(gene_id_style) && is.null(cur$microarray)) {
    gene_id_style <- env$gene_id
  }
  # Take the first file in the directory by default
  dir <- here(env$raw_path, name)
  if (!is.null(lget(cur, "files", list())$raw)) {
    file <- here(dir, cur$files$raw)
  } else {
    file <- list.files(dir, full.names = TRUE)[1]
  }

  expr <- suppressMessages(
    if (str_ends(file, "csv")) read_csv(file) else read_tsv(file)
  )
  meta <- dplyr::filter(mdata, dataset == name & join_id %in% colnames(expr))
  expr <- distinct(expr, across(1), .keep_all = TRUE)
  expr <- column_to_rownames(expr, var = colnames(expr)[1])
  colnames(expr) <- paste0(name, "_", colnames(expr))
  if (!is.null(cur$microarray)) {
    platform <- unique(meta$platform)
    # Will collapse probes mapping to the same genes with the mean, if any exist,
    # following the protocol in the I-SPY 2 dataset (GSE194040)
    if (length(platform) != 1) {
      stop("There should only be one platform per dataset")
    }
    gpl <- getGEO(platform, destdir = geo_cache) |> Table()
    print(platform)
    print(head(gpl))
    id_col <- if (is.null(cur$probeset_id_col)) "ID" else cur$probeset_id_col
    probeid2symbol <- setNames(gpl$GeneName, gpl[[id_col]])
    to_symbol <- probeid2symbol[rownames(expr)]
    if (any(duplicated(to_symbol))) {
      expr <- as_tibble(expr) |>
        mutate(ID = to_symbol) |>
        filter(!is.na(ID)) |>
        group_by(ID) |>
        summarise(across(everything(), mean)) |>
        column_to_rownames(var = "ID")
    } else {
      mask <- !is.na(to_symbol)
      expr <- expr[mas, ]
      rownames(expr) <- to_symbol[mask]
    }
    gene_id_style <- "SYMBOL"
  }

  if (gene_id_style != "ncbi") {
    kt <- gene_id_style
    mapped <- mapIds(db, rownames(expr), keytype = kt, column = "ENTREZID") |>
      unlist()
    mask <- !duplicated(mapped) & !is.na(mapped)
    expr <- expr[mask, ]
    rownames(expr) <- mapped[mask]
  }

  meta <- AnnotatedDataFrame(column_to_rownames(meta, var = "join_id"))
  rownames(meta) <- paste0(name, "_", rownames(meta))
  ## message(glue("n samples in meta {nrow(meta)}"))
  ## message(glue("n samples in expr {ncol(expr)}"))
  expr <- expr[, colnames(expr) %in% rownames(meta)]
  meta <- meta[sort(rownames(meta)), ] # Required to construct eset
  expr <- expr[, sort(colnames(expr))] |> as.matrix()
  Biobase::ExpressionSet(assayData = expr, phenoData = meta)
}) |>
  `names<-`(names(env$datasets))

filtered_mdata <- lapply(esets, \(x) as_tibble(pData(x))) |> bind_rows()

write_tsv(filtered_mdata, here(outdir, glue("sample_metadata-{date}.tsv")))
