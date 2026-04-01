suppressMessages({
  library(shiny)
  library(optparse)
  library(bslib)
  library(tidyverse)
  library(polars)
  library(checkmate)
  library(gt)
  library(glue)
})

read_gs <- function(file) {
  if (str_ends(file, "yaml")) {
    lst <- yaml::read_yaml(file)
    pl$LazyFrame(gs = names(lst), gene = lst)$explode("gene")
  } else {
    if (str_ends(file, "csv")) {
      sep <- ","
    } else {
      sep <- "\t"
    }
    pl$scan_csv(file, separator = sep, has_header = FALSE)$rename(
      column_1 = "gs",
      column_2 = "gene"
    )
  }
}


## * Table functions

make_de_table <- function(lfs, gene_col, gene_group_col = NULL, cfg, input) {
  lf <- lfs[[input$gs_name]]$filter(pl$col("is_de"))$with_columns(pl$col(
    "gs"
  )$list$join(", "))$drop("is_de")
  if (!is.null(gene_group_col) && !is.null(input$gene_group)) {
    df <- lf$filter(pl$col(gene_group_col) == input$gene_group)
  }
  lf <- lf$filter(pl$col(cfg$lfc_column)$abs() > input$lfc_thresh)
  if (!is.null(cfg$significance_column)) {
    lf <- lf$filter(pl$col(cfg$significance_column) > input$sig_thresh)
  }

  if (input$lfc_direction == "Positive") {
    lf <- lf$filter(pl$col(cfg$lfc_column) > 0)
  } else if (input$lfc_direction == "Negative") {
    lf <- lf$filter(pl$col(cfg$lfc_column) < 0)
  }
  gt(as_tibble(lf)) |>
    fmt_number(columns = cfg$other_numeric) |>
    opt_interactive(page_size_default = 10, use_search = TRUE)
}

make_gs_table <- function(lfs, cfg, input) {
  lfc_col <- cfg$lfc_column
  lf <- lfs[[input$gs_name]]
  grouped <- lf$explode("gs")$group_by("gs")$agg(
    pl$len()$alias("Size"),
    pl$col("is_de")$sum()$alias("n DE"),
    pl$col(lfc_col)$max()$alias("Max LFC"),
    pl$col(lfc_col)$mean()$alias("Mean LFC")
  )$with_columns((pl$col("n DE") / pl$col("Size"))$alias("% DE"))
  numeric_cols <- names(grouped)[-1]
  gt(as_tibble(grouped)) |>
    fmt_number(columns = numeric_cols) |>
    opt_interactive(page_size_default = 10, use_search = TRUE)
}

set_up_cfg <- function() {
  parser <- OptionParser()
  parser <- add_option(
    parser,
    c("-i", "--input"),
    type = "character",
    help = "Input filename"
  )
  parser <- add_option(
    parser,
    c("-c", "--config"),
    type = "character",
    help = "Configuration file",
    default = NULL
  )
  parser <- add_option(
    parser,
    c("-p", "--port"),
    type = "character",
    help = "Shiny port to run server on"
  )
  args <- parse_args(parser)
  if (!is.null(args$config)) {
    cfg <- yaml::read_yaml(args$config)
    assert_named(cfg)
    assert_names(
      names(cfg),
      subset.of = c(
        "input",
        "gene_sets",
        "gene_column",
        "gene_name_format",
        "grouping_column",
        "lfc_column",
        "significance_column",
        "other_numeric"
      )
    )
  } else {
    cfg <- list(
      gene_name_format = "symbol",
      lfc_column = "lfc",
      significance_column = "p_value"
    )
  }
  assert_choice(cfg$gene_name_format, c("symbol", "ensembl", "ncbi"))
  assert_list(cfg$gene_sets)
  assert_named(cfg$gene_sets)
  list(args = args, cfg = cfg)
}


## * Set up input files

tmp <- set_up_cfg()
args <- tmp$args
cfg <- tmp$cfg

input <- args$input %||% cfg$input
g_name_format <- cfg$gene_name_format
g_col <- cfg$gene_column

sep <- ifelse(str_ends(input, "csv"), ",", "\t")

cols_select <- c(
  g_col,
  cfg$lfc_column,
  cfg$significance_column,
  cfg$grouping_column,
  cfg$other_numeric
)

de <- pl$scan_csv(input, separator = sep)$select(cols_select)

cols_order <- c(cols_select, "gs", "is_de")

all_de_genes <- list(de$collect()[[g_col]]$unique()$to_r_vector())

get_gs <- function(x) {
  lf <- read_gs(x)
  gs_keep <- lf$filter(pl$col("gene")$is_in(all_de_genes))$collect()[[
    "gs"
  ]]$unique()$to_r_vector()
  gs_keep <- list(gs_keep)

  lf$filter(pl$col(
    "gs"
  )$is_in(gs_keep))$group_by("gene")$agg(pl$col("gs"))$with_columns(pl$col(
    "gene"
  )$is_in(all_de_genes)$alias(
    "is_de"
  ))$join(
    de,
    how = "full",
    right_on = "gene",
    left_on = g_col
  )$select(cols_order)
}

with_gs <- lapply(cfg$gene_sets, get_gs)

gs_names <- names(with_gs)

g_group_col <- cfg$grouping_column
if (!is.null(g_group_col)) {
  g_groups <- de$collect()[[g_group_col]]$unique()$to_r_vector()
} else {
  g_groups <- character(0)
}


## * Application

## ** UI

ui <- page_navbar(
  id = "nav",
  sidebar = sidebar(
    selectInput(
      inputId = "gs_name",
      label = "Gene set definitions",
      choices = gs_names
    ),
    conditionalPanel(
      condition = !is.null(g_group_col),
      selectInput(
        inputId = "gene_group",
        label = "Group",
        choices = g_groups
      )
    ),
    numericInput(
      inputId = "lfc_thresh",
      label = "LFC threshold (absolute)",
      value = 0
    ),
    selectInput(
      inputId = "lfc_direction",
      label = "LFC direction",
      choices = c("Positive", "Negative", "Both")
    ),
    conditionalPanel(
      condition = !is.null(cfg$significance_column),
      numericInput(
        inputId = "sig_thresh",
        label = "Significance threshold",
        value = 0.05
      )
    )
  ),
  nav_panel(
    "DE genes",
    gt_output(outputId = "de_table"),
  ),
  nav_panel(
    "Gene set statistics",
    gt_output(outputId = "gs_table"),
  )
)


## ** Server
server <- function(input, output, session) {
  output$de_table <- render_gt(
    expr = make_de_table(
      lfs = with_gs,
      gene_col = g_col,
      gene_group_col = g_group_col,
      input = input,
      cfg = cfg
    )
  )
  output$gs_table <- render_gt(
    expr = make_gs_table(lfs = with_gs, cfg = cfg, input = input)
  )
}


app <- shinyApp(ui = ui, server = server)
runApp(app, launch.browser = TRUE, port = 4214)
