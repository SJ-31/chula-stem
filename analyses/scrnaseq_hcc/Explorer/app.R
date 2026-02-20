suppressMessages({
  library(shiny)
  library(bslib)
  library(reticulate)
  library(tidyverse)
  library(polars)
  library(tidypolars)
  library(gt)
  library(memoise)
  library(glue)
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

# TODO: you can try to modularize this

## * Set up data
anno_dir <- glue("{env$outdir}/annotations")


scan_if_exists <- function(file) {
  if (file.exists(file)) {
    pl$scan_csv(file)
  } else {
    NULL
  }
}

get_group_map <- function(lf, key_col, val_col) {
  c <- lf$unique(c(key_col, val_col))$group_by(key_col)$agg(
    !!!list(pl$col(val_col))
  )$collect()
  as.list(c[[val_col]]) |> `names<-`(as.character(c[[key_col]]))
}


cc_de_edger <- scan_if_exists(glue("{anno_dir}/clusters-edgeR_de.csv"))

if (!is.null(cc_de_edger)) {
  clusterings2contrasts <- get_group_map(cc_de_edger, "clustering", "contrast")
} else {
  clusterings2contrasts <- NULL
}

cc_markers <- scan_if_exists(glue("{anno_dir}/marker_gene_activity.csv"))$cast(
  group = pl$String
)
cc_gs <- scan_if_exists(glue("{anno_dir}/gene_set_activity.csv"))$cast(
  group = pl$String
)
clusterings2cluster_names <- get_group_map(cc_markers, "clustering", "group")

cc_clusterings <- NULL

## * Table functions

## ** formatting

# TODO: make a unified format between scVI and edgeR

## ** GT

make_table <- function(
  ltb,
  clustering_name,
  group,
  pos,
  group_col = "group",
  title,
  direction_col = "stat",
  number_cols = c("stat", "meanchange", "pval", "padj"),
  cols_remove = c("group", "clustering", "reference")
) {
  ltb$filter(
    pl$col(group_col) == group,
    pl$col("clustering") == clustering_name
  ) |>
    arrange(desc(abs(!!as.symbol(direction_col)))) |>
    filter(
      (pos & !!as.symbol(direction_col) > 0) |
        (!pos & !!as.symbol(direction_col) < 0)
    ) |>
    select(-any_of(cols_remove)) |>
    collect() |>
    gt() |>
    fmt_number(columns = number_cols) |>
    opt_interactive(page_size_default = 5, use_search = TRUE) |>
    tab_header(title)
}


## make_table <- (make_table_helper)

## * Application

provide_clustering_choice <- function(lf) {
  if (is.null(lf)) {
    choices <- character(0)
  } else {
    if (is.null(cc_clusterings)) {
      cc_clusterings <<- lf |>
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
    gt_output(outputId = "cc_enrich_gs_table"),
    gt_output(outputId = "cc_enrich_markers_table")
  ),
  nav_panel(
    "Cell cluster DE",
    gt_output(outputId = "cc_de_edger_table"),
    gt_output(outputId = "cc_de_scvi_table"),
    gt_output(outputId = "cc_de_scvi_extra_table")
  ),
  nav_panel("Sample-level DE")
)

## * Server
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
  ## ** CC enrichment panel
  output$cc_enrich_gs_table <- render_gt(
    expr = make_table(
      cc_gs,
      input$clustering_choice,
      group = input$cluster_name,
      pos = input$stat_positive,
      title = "Pathways"
    )
  )
  output$cc_enrich_markers_table <- render_gt(
    expr = make_table(
      cc_markers,
      input$clustering_choice,
      group = input$cluster_name,
      pos = input$stat_positive,
      title = "Cell markers"
    )
  )
  ## ** CE DE panel
  output$cc_de_edger_table <- render_gt(
    expr = make_table(
      cc_de_edger,
      input$clustering_choice,
      input$contrast,
      input$stat_positive,
      group_col = "contrast",
      direction_col = "logFC",
      number_cols = c("logFC", "unshrunk.logFC", "logCPM", "PValue", "FDR"),
      cols_remove = c("clustering", "contrast")
    )
  )
}

shinyApp(ui = ui, server = server)
