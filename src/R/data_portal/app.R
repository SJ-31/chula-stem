suppressMessages({
  options(shiny.autoreload = TRUE)
  library(shiny)
  library(bslib)
  library(tidyverse)
  library(googlesheets4)
  library(reactable)
})

link <- "https://docs.google.com/spreadsheets/d/1h8aGAyhQk1IL2MALBZqmARUjvurcYIdCfi2GcHPQVc4/edit?gid=870277624#gid=870277624"


setup_data <- function() {
  data <- list(cohorts = list(), ttypes = list())
  df <- read_sheet(link) |>
    mutate(tumor_type = str_to_upper(tumor_type)) |>
    select(-date_received) |>
    relocate(path, .after = everything())
  modalities <- as.character(df[["modality"]]) |> unique()
  data$dfs <- lapply(modalities, \(x) {
    filtered <- df |> filter(modality == x)
    data$selections[[x]] <<- select(df, c("cohort", "tumor_type")) |> distinct()
    filtered |> select(-modality)
  }) |>
    `names<-`(modalities)
  data
}

D <- setup_data()


clean_modality <- function(name) {
  str_replace_all(str_to_lower(name), "-", "_")
}

get_selection <- function(input, condition_on, out_col) {
  choice_df <- D$selections[[clean_modality(input$nav)]]

  if (input[[condition_on]] != "-") {
    mask <- choice_df[[condition_on]] == input[[condition_on]]
    filtered <- choice_df[mask, ]
  } else {
    filtered <- choice_df
  }
  sort(c(unique(as.character(filtered[[out_col]])), "-"))
}

rename_df <- function(df) {
  df <- df |>
    rename(
      Case = "case_name",
      Path = "path",
      PBMC = "has_pbmc",
      Tumor = "has_tumor",
      Raw = "has_raw",
      Processed = "has_processed"
    )
  if ("tumor_type" %in% colnames(df)) {
    df <- df |> rename("Tumor Type" = "tumor_type")
  }
  if ("cohort" %in% colnames(df)) {
    df <- df |> rename(Cohort = "cohort")
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

  cohort_val <- input$cohort
  ttype_val <- input$tumor_type

  if (cohort_val != "-" && ttype_val != "-") {
    filtered <- df |>
      filter(
        cohort == cohort_val & tumor_type == ttype_val
      ) |>
      select(-cohort, -tumor_type)
  } else if (cohort_val != "-") {
    filtered <- df |>
      filter(cohort == cohort_val) |>
      select(-cohort)
  } else if (ttype_val != "-") {
    filtered <- df |>
      filter(tumor_type == ttype_val) |>
      select(-tumor_type)
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
