suppressMessages({
  options(shiny.autoreload = TRUE)
  library(shiny)
  library(bslib)
  library(tidyverse)
  library(googlesheets4)
  library(reactable)
  library(crosstalk)
})
gs4_deauth()

source("setup.R")

ui <- page_navbar(
  theme = bs_theme(bootswatch = "flatly"),
  id = "nav",
  nav_panel("Clinical", reactableOutput("clin_tab")),
  nav_panel("Metadata", reactableOutput("meta_tab")),
  nav_panel(
    "Catalog",
    reactableOutput("mtab")
  ),
  sidebar = sidebar(
    h3("Filters"),
    filter_checkbox("modality", "Modality", D$all, ~Modality),
    filter_select("cohort", "Cohort", D$all, ~Cohort),
    filter_checkbox("ttype", "Tumor Type", D$all, ~`Tumor Type`),
    csv_download_button("mtab")
  )
)

server <- function(input, output, session) {
  output$mtab <- renderReactable(make_sample_table())
  output$clin_tab <- renderReactable(make_clinical_table())
  output$meta_tab <- renderReactable(make_metadata_table())
  selected <- reactive(getReactableState("mtab", "selected"))
}

shinyApp(ui = ui, server = server)
