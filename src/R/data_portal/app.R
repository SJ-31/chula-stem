suppressMessages({
  options(shiny.autoreload = TRUE)
  library(shiny)
  library(bslib)
  library(tidyverse)
  library(polars)
  library(reactable)
  library(glue)
})

sample_file <- "/home/shannc/Downloads/2026-05-07-record_cases-autogen.csv"


setup_data <- function(manifest) {
  data <- list(cohorts = list(), ttypes = list())
  df <- pl$read_csv(manifest)
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
  c(unique(as.character(filtered[[out_col]])), "-")
}

make_sample_table <- function(input) {
  df <- D$dfs[[clean_modality(input$nav)]]

  print(df)

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
  reactable(as.data.frame(filtered))
}

## plot_samples <- function(df) {

## }

ui <- page_navbar(
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
  )
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
