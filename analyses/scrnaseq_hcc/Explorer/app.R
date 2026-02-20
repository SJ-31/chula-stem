suppressMessages({
  library(shiny)
  library(bslib)
  library(reticulate)
  library(tidyverse)
  library(duckdb)
  library(gt)
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


cc_de <- NULL
clusterings2contrasts <- list()
try({
  cc_de <- tbl(con, glue("{anno_dir}/clusters-edgeR_de.csv"))
  clusterings2contrasts <- cc_de |>
    group_by(clustering) |>
    collect() |>
    summarise(name = list(unique(contrast))) |>
    deframe()
})


cc_markers <- tbl(con, glue("read_csv('{anno_dir}/marker_gene_activity.csv')"))
cc_gs <- tbl(con, glue("{anno_dir}/gene_set_activity.csv"))
clusterings2cluster_names <- cc_gs |>
  group_by(clustering) |>
  collect() |>
  summarise(name = list(unique(group))) |>
  deframe()


cc_clusterings <- NULL

## * Table functions

make_enrichment_table <- function(
  ltb,
  clustering,
  cluster,
  pos,
  title,
  direction_col = "stat"
) {
  ltb |>
    filter(clustering == clustering & group == cluster) |>
    filter((pos & stat > 0) | (!pos & stat < 0)) |>
    select(-group, -clustering, -reference) |>
    collect() |>
    gt() |>
    fmt_number(columns = c(stat, meanchange, pval, padj)) |>
    opt_interactive() |>
    tab_header(title)
}


## * Application

provide_clustering_choice <- function(lazy_tb) {
  if (is.null(lazy_tb)) {
    choices <- character(0)
  } else {
    if (is.null(cc_clusterings)) {
      cc_clusterings <<- lazy_tb |>
        distinct(clustering) |>
        collect() |>
        pluck("clustering")
    }
    choices <- cc_clusterings
  }
  selectInput(
    inputId = "clustering_choice",
    label = "Clustering",
    choices = choices
  )
}


ui <- page_navbar(
  id = "nav",
  sidebar = sidebar(
    conditionalPanel(
      condition = "input.nav === 'Cell cluster enrichment'",
      provide_clustering_choice(cc_gs),
      selectInput(
        inputId = "cluster_name",
        label = "Cluster",
        choices = character(0)
      ),
      checkboxInput("stat_positive", "Enriched relative to others", TRUE),
      "The results show groups enriched in the chosen cluster against all other clusters (OvR)"
    ),
    conditionalPanel(
      condition = "input.nav === 'Cell cluster DE'",
      provide_clustering_choice(cc_gs),
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
  nav_panel(
    "Cell cluster enrichment",
    gt_output(outputId = "cc_enrich_gs"),
    gt_output(outputId = "cc_enrich_markers")
  ),
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
  output$cc_enrich_gs <- render_gt(
    expr = make_enrichment_table(
      cc_gs,
      input$clustering_choice,
      input$cluster_name,
      input$stat_positive,
      "Pathways"
    )
  )
  output$cc_enrich_markers <- render_gt(
    expr = make_enrichment_table(
      cc_markers,
      input$clustering_choice,
      input$cluster_name,
      input$stat_positive,
      "Cell markers"
    )
  )
}

shinyApp(ui = ui, server = server)
