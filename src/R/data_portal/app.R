suppressMessages({
  options(shiny.autoreload = TRUE)
  library(shiny)
  library(bslib)
  library(tidyverse)
  library(polars)
  library(reactable)
})

link <- "https://docs.google.com/spreadsheets/d/1h8aGAyhQk1IL2MALBZqmARUjvurcYIdCfi2GcHPQVc4/edit?gid=870277624#gid=870277624"

read_sheet(
  link,
  col_types = "c"
)

sample_file <- "/home/shannc/Downloads/2026-05-07-record_cases-autogen.csv"


setup_data <- function(manifest) {
  data <- list(cohorts = list(), ttypes = list())
  df <- pl$read_csv(manifest)$with_columns(pl$col(
    "tumor_type"
  )$str$to_uppercase())$drop("date_received")
  cols <- c(discard(df$columns, \(x) x == "path"), "path")
  df <- df$select(cols)
  modalities <- as.character(df[["modality"]]) |> unique()
  data$dfs <- lapply(modalities, \(x) {
    filtered <- df$filter(pl$col("modality") == x)
    data$selections[[x]] <<- filtered$select(c("cohort", "tumor_type"))$unique()
    filtered$drop("modality")
  }) |>
    `names<-`(modalities)
  data
}

D <- setup_data(sample_file)


clean_modality <- function(name) {
  str_replace_all(str_to_lower(name), "-", "_")
}

get_selection <- function(input, condition_on, out_col) {
  choice_df <- D$selections[[clean_modality(input$nav)]]

  if (input[[condition_on]] != "-") {
    filtered <- choice_df$filter(
      pl$col(condition_on) == input[[condition_on]]
    )
  } else {
    filtered <- choice_df
  }
  sort(c(unique(as.character(filtered[[out_col]])), "-"))
}

rename_df <- function(df) {
  df <- df$rename(
    case_name = "Case",
    path = "Path",
    has_pbmc = "PBMC",
    has_tumor = "Tumor",
    has_raw = "Raw",
    has_processed = "Processed"
  )
  if ("tumor_type" %in% df$columns) {
    df <- df$rename(tumor_type = "Tumor Type")
  }
  if ("cohort" %in% df$columns) {
    df <- df$rename(cohort = "Cohort")
  }
  df
}

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

make_sample_table <- function(input) {
  df <- D$dfs[[clean_modality(input$nav)]]

  cohort <- input$cohort
  tumor_type <- input$tumor_type

  if (cohort != "-" && tumor_type != "-") {
    filtered <- df$filter(
      pl$col("cohort") == cohort & pl$col("tumor_type") == tumor_type
    )$drop(c("cohort", "tumor_type"))
  } else if (cohort != "-") {
    filtered <- df$filter(pl$col("cohort") == cohort)$drop("cohort")
  } else if (tumor_type != "-") {
    filtered <- df$filter(pl$col("tumor_type") == tumor_type)$drop("tumor_type")
  } else {
    filtered <- df
  }
  reactable(
    as.data.frame(rename_df(filtered)),
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

## plot_samples <- function(df) {

## }

ui <- page_navbar(
  theme = bs_theme(bootswatch = "flatly"),
  id = "nav",
  nav_panel("Exome", reactableOutput("t1")),
  nav_panel("RNA-seq", reactableOutput("t2")),
  nav_panel("scRNA-seq", reactableOutput("t3")),
  nav_panel("TCR-seq", reactableOutput("t4")),
  nav_panel("sc-ATAC-seq", reactableOutput("t5")),
  sidebar = sidebar(
    selectInput(
      inputId = "cohort",
      selected = "-",
      label = "Cohort",
      choices = "-"
    ),
    selectInput(
      inputId = "tumor_type",
      selected = "-",
      label = "Tumor Type",
      choices = "-"
    ),
  ),
  title = "Modalities"
)

server <- function(input, output, session) {
  output$t1 <- renderReactable(make_sample_table(input))
  output$t2 <- renderReactable(make_sample_table(input))
  output$t3 <- renderReactable(make_sample_table(input))
  output$t4 <- renderReactable(make_sample_table(input))
  output$t5 <- renderReactable(make_sample_table(input))
  observeEvent(input$nav, {
    updateSelectInput(
      inputId = "cohort",
      selected = "-",
      choices = get_selection(input, "tumor_type", "cohort")
    )
    updateSelectInput(
      inputId = "tumor_type",
      selected = "-",
      choices = get_selection(input, "cohort", "tumor_type")
    )
  })
  observeEvent(input$cohort, {
    choices <- get_selection(input, "cohort", "tumor_type")
    if (input$tumor_type %notin% choices) {
      selected <- "-"
    } else {
      selected <- input$tumor_type
    }
    updateSelectInput(
      inputId = "tumor_type",
      selected = selected,
      choices = choices
    )
  })
  observeEvent(input$tumor_type, {
    choices <- get_selection(input, "tumor_type", "cohort")
    if (input$cohort %notin% choices) {
      selected <- "-"
    } else {
      selected <- input$cohort
    }
    updateSelectInput(
      inputId = "cohort",
      selected = selected,
      choices = get_selection(input, "tumor_type", "cohort")
    )
  })
}

app <- shinyApp(ui = ui, server = server)
runApp(app, launch.browser = TRUE, port = 4214)
