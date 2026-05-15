library(tidyverse)
library(crosstalk)
library(checkmate)

##  Validate config
local({
  data_sheets <- config::get("data_sheets")
  assert_subset(
    names(data_sheets),
    choices = c("samples", "other")
  )
  assert_string(data_sheets$samples)
  assert_string(data_sheets$other)
})


JOIN_ON <- c(
  "tumor_type",
  "case_name",
  "cohort"
)

EXPECTED_COLS <- list(
  clinical = c("HN", "diagnosis", "note"),
  metadata = c(
    "platform",
    "ngs_provider",
    "sample_collection_date",
    "received_date"
  )
)

read_from_other <- function(link, sheet_name, grouped_sample_df) {
  read_sheet(link, sheet = sheet_name) |>
    mutate(tumor_type = str_to_upper(tumor_type)) |>
    left_join(grouped_sample_df, by = JOIN_ON) |>
    mutate(`Has data?` = ifelse(is.na(Modality), "F", "T")) |>
    select(-Modality)
}


rename_df <- function(df) {
  lookup <- c(
    Case = "case_name",
    Path = "path",
    PBMC = "has_pbmc",
    Tumor = "has_tumor",
    Raw = "has_raw",
    Processed = "has_processed",
    Modality = "modality",
    Diagnosis = "diagnosis",
    `NGS platform` = "platform",
    `NGS provider` = "ngs_provider",
    `Sample collection` = "sample_collection_date",
    `Received data` = "received_date",
    `Tumor Type` = "tumor_type",
    Cohort = "cohort",
    Note = "note"
  )
  df |> rename(any_of(lookup))
}

D <- local({
  data <- list()
  sheets <- config::get("data_sheets")
  other_sheet <- sheets$other
  df <- read_sheet(sheets$samples) |>
    mutate(
      tumor_type = str_to_upper(tumor_type),
      modality = replace_values(
        modality,
        "exome" ~ "Exome",
        "rna_seq" ~ "RNA seq",
        "scrna_seq" ~ "scRNA seq",
        "sc_atac_seq" ~ "scATAC seq",
        "tcr_seq" ~ "TCR seq"
      ),
    ) |>
    select(-date_received) |>
    relocate(path, .after = everything())
  grouped <- df |>
    group_by(cohort, case_name, tumor_type) |>
    summarise(Modality = paste0(modality, collapse = "; "))
  data$all <- rename_df(df) |> SharedData$new(group = "tables")
  data$clinical <- read_from_other(other_sheet, "clinical", grouped) |>
    rename_df() |>
    SharedData$new(group = "tables")
  data$meta <- read_from_other(other_sheet, "metadata", grouped) |>
    rename_df() |>
    SharedData$new(group = "tables")
  data
})


select_filter <- function(values, name) {
  tags$select(
    onchange = sprintf(
      "Reactable.setFilter('main_tab', '%s', event.target.value || undefined)",
      name
    ),
    tags$option(value = "", "All"),
    lapply(unique(values), tags$option)
  )
}

bool_col <- colDef(
  filterable = TRUE,
  filterInput = select_filter,
  align = "center",
  cell = \(v) {
    if (v == "F") "\u274c" else "\u2714\ufe0f"
  }
)

csv_download_button <- function(
  id,
  filename = "data.csv",
  label = "Download as CSV"
) {
  tags$button(
    tagList(icon("download"), label),
    onclick = sprintf("Reactable.downloadDataCSV('%s', '%s')", id, filename)
  )
}

make_clinical_table <- function(input) {
  reactable(
    D$clinical,
    searchable = TRUE,
    wrap = FALSE,
    columns = list(
      Note = colDef(
        cell = function(value, index, name) "",
        html = TRUE,
        details = JS(
          "function(rowInfo) {
        return `</br>${rowInfo.values['Note']}</br>`
}"
        ),
      ),
      `Has data?` = bool_col
    )
  )
}

make_metadata_table <- function() {
  reactable(
    D$meta,
    columnGroups = list(
      colGroup(
        name = "Dates",
        columns = c("Sample collection", "Received data")
      )
    ),
    columns = list(
      `Sample collection` = colDef(format = colFormat(date = TRUE)),
      `Received data` = colDef(format = colFormat(date = TRUE))
    )
  )
}


make_sample_table <- function() {
  reactable(
    D$all,
    pagination = FALSE,
    searchable = TRUE,
    showPagination = TRUE,
    resizable = TRUE,
    wrap = FALSE,
    columns = list(
      PBMC = bool_col,
      Raw = bool_col,
      Processed = bool_col,
      Tumor = bool_col,
      Path = colDef(filterable = FALSE)
    ),
    elementId = "main_tab",
    columnGroups = list(
      colGroup(
        name = "Available Data",
        columns = c("Processed", "Raw")
      ),
      colGroup(
        name = "Available Samples",
        columns = c("Tumor", "PBMC")
      )
    )
  )
}
