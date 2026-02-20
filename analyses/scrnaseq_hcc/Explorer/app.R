suppressMessages({
  library(shiny)
  library(bslib)
  library(reticulate)
  library(tidyverse)
  library(duckdb)
  library(glue)
  library(dbplyr)
  use_condaenv("stem-base")
})
library(here)

yte <- import("yte")
env <- yte$process_yaml(paste0(
  read_lines(
    "/home/shannc/Bio_SDD/chula-stem/analyses/scrnaseq_hcc/cellranger_config.yaml"
  ),
  # TODO: make this a cli parameter
  collapse = "\n"
))

## * Set up data
anno_dir <- glue("{env$outdir}/annotations")

con <- dbConnect(duckdb())

cc_gs_markers <- tbl(con, glue("read_csv('{anno_dir}/*_activity.csv')"))
cc_de <- tbl(con, glue("{anno_dir}/clusters-edgeR_de.csv"))

clusterings2cluster_names <- cc_gs_markers |>
  group_by(clustering) |>
  collect() |>
  summarise(name = list(unique(group))) |>
  deframe()
clusterings2contrasts <- cc_de |>
  group_by(clustering) |>
  collect() |>
  summarise(name = list(unique(contrast))) |>
  deframe()

cc_clusterings <- NULL

## * Application

provide_clustering_choice <- function(lazy_tb) {
  if (is.null(cc_clusterings)) {
    cc_clusterings <<- lazy_tb |>
      distinct(clustering) |>
      collect() |>
      pluck("clustering")
  }
  selectInput(
    inputId = "clustering_choice",
    label = "Clustering",
    choices = cc_clusterings
  )
}


ui <- page_navbar(
  id = "nav",
  sidebar = sidebar(
    conditionalPanel(
      condition = "input.nav === 'Cell cluster enrichment'",
      provide_clustering_choice(cc_gs_markers),
      selectInput(
        inputId = "cluster_name",
        label = "Cluster",
        choices = character(0)
      ),
      "The results show groups enriched in the chosen cluster against all other clusters (OvR)"
    ),
    conditionalPanel(
      condition = "input.nav === 'Cell cluster DE'",
      provide_clustering_choice(cc_gs_markers),
      selectInput(
        inputId = "contrast",
        label = "Contrast",
        choices = character(0)
      )
    ),
    conditionalPanel(
      condition = "input.nav === 'Sample-level DE'"
    )
  ),
  nav_panel("Cell cluster enrichment", "FOO"),
  nav_panel("Cell cluster DE"),
  nav_panel("Sample-level DE")
)

server <- function(input, output, session) {
  observeEvent(input$clustering_choice, {
    if (input$nav == "Cell cluster enrichment") {
      updateSelectInput(
        session,
        "cluster_name",
        choices = clusterings2cluster_names[[input$clustering_choice]]
      )
    } else {
      updateSelectInput(
        session,
        "contrast",
        choices = clusterings2contrasts[[input$clustering_choice]]
      )
    }
  })
}

shinyApp(ui = ui, server = server)
