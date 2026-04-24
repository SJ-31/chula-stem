suppressMessages({
  library(shiny)
  library(bslib)
  library(reticulate)
  library(tidyverse)
  library(polars)
  library(checkmate)
  library(tidypolars)
  library(gt)
  library(memoise)
  library(shinyjqui)
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


cc_de_edger <- scan_if_exists(glue("{anno_dir}/clusters-edgeR_de.csv"))$cast(
  contrast = pl$String
)
cc_de_scvi <- scan_if_exists(glue("{anno_dir}/clusters-scVI_de.csv"))$cast(
  contast = pl$String
)

if (!is.null(cc_de_edger)) {
  clusterings2contrasts <- get_group_map(cc_de_edger, "clustering", "contrast")
} else {
  clusterings2contrasts <- NULL
}

edger_cols <- c(
  "logFC",
  "unshrunk.logFC",
  "logCPM",
  "PValue",
  "FDR",
  "lfc"
)

cc_markers <- scan_if_exists(glue("{anno_dir}/marker_gene_activity.csv"))$cast(
  group = pl$String
)$filter(!pl$col("redundant"))
cc_gs <- scan_if_exists(glue("{anno_dir}/gene_set_activity.csv"))$cast(
  group = pl$String
)$filter(!pl$col("redundant"))

sample_de_tmp <- scan_if_exists(glue("{anno_dir}/samples_de.csv"))
sample_de_edger <- sample_de_tmp$filter(pl$col(
  "analysis_group"
)$str$starts_with(
  "edgeR"
))$select(edger_cols)
sample_de_scvi <- sample_de_tmp$filter(pl$col(
  "analysis_group"
)$str$starts_with(
  "scVI"
))$drop(edger_cols)

sample_enrich <- scan_if_exists(glue("{anno_dir}/samples_de_gprofiler.csv"))

clusterings2cluster_names <- get_group_map(cc_markers, "clustering", "group")

cc_clusterings <- NULL
sample_analysis_groups <- NULL

## * Table functions

## ** Formatting

# TODO: make a unified format between scVI and edgeR

## ** GT

make_table <- function(
  ltb,
  subgroup_name,
  group,
  pos,
  group_col = "group",
  subgroup_col = "clustering",
  title,
  direction_col = "stat",
  number_cols = c("stat", "meanchange", "pval", "padj"),
  cols_remove = c("group", "clustering", "reference")
) {
  ltb$filter(
    pl$col(group_col) == group,
    pl$col(subgroup_col) == subgroup_name
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

provide_group_choice <- function(lf, for_clusters = TRUE) {
  if (is.null(lf)) {
    choices <- character(0)
  } else if (for_clusters) {
    var_name <- "cc_clusterings"
    colname <- "clustering"
    input_id <- "clustering_choice"
    label <- "Clustering"
  } else {
    var_name <- "sample_analysis_groups"
    colname <- "analysis_group"
    input_id <- "agroup_choice"
    label <- "Analysis Group"
  }
  var_vals <- globalenv()[[var_name]]
  if (is.null(var_vals)) {
    assign(
      var_name,
      lf |>
        distinct(!!as.symbol(colname)) |>
        collect() |>
        pluck(colname),
      envir = globalenv()
    )
  }
  choices <- globalenv()[[var_name]]
  selectInput(
    inputId = input_id,
    label = label,
    choices = choices
  )
}

ui <- page_navbar(
  id = "nav",
  sidebar = sidebar(
    conditionalPanel(
      condition = "input.nav === 'Cell cluster enrichment'",
      provide_group_choice(cc_gs),
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
      provide_group_choice(cc_gs),
      selectInput(
        inputId = "contrast",
        label = "Contrast",
        choices = character(0)
      )
    ),
    conditionalPanel(
      condition = "input.nav === 'Sample-level DE'",
      provide_group_choice(sample_de_scvi, FALSE),
      selectInput(
        inputId = "contrast",
        label = "Contrast",
        choices = character(0)
      )
    ),
    conditionalPanel(
      condition = "input.nav === 'Sample-level enrichment'",
      provide_group_choice(sample_de_scvi, FALSE),
      selectInput(
        inputId = "contrast",
        label = "Contrast",
        choices = character(0)
      )
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
  nav_panel(
    "Sample-level DE",
    gt_output(outputId = "sample_de_edger_table"),
    gt_output(outputId = "sample_de_scvi_table")
  ),
  nav_panel(
    "Sample-level enrichment",
    gt_output(outputId = "sample_enrich_table")
  )
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
      subgroup_name = input$clustering_choice,
      group = input$contrast,
      pos = input$stat_positive,
      group_col = "contrast",
      direction_col = "logFC",
      number_cols = c("logFC", "unshrunk.logFC", "logCPM", "PValue", "FDR"),
      cols_remove = c("clustering", "contrast"),
      title = "edgeR"
    )
  )
  output$cc_de_scvi_table <- render_gt(
    expr = make_table(
      cc_de_scvi,
      subgroup_name = input$clustering_choice,
      group = input$contrast,
      pos = input$stat_positive,
      group_col = "contrast",
      direction_col = "logFC",
      number_cols = c("proba_de"),
      cols_remove = c(
        "clustering",
        "contrast",
        "proba_not_de",
        "group1",
        "group2"
      ),
      title = "scVI",
    ) |>
      fmt_number(starts_with("lfc"), starts_with("raw"))
  )
  ## TODO: add these in [2026-03-02 Mon]
  output$sample_de_edger_table <- render_gt(
    expr = make_table(
      sample_de_edger,
      subgroup_name = input$comparison_choice,
      group = input$contrast,
      pos = input$stat_positive,
      subgroup_col = "analysis_group",
      direction_col = "logFC",
      title = "edgeR"
    )
  )
  output$sample_de_scvi_table <- render_gt(expr = make_table())
  ## output$sample_enrich_table <- make_table()
}

shinyApp(ui = ui, server = server)
