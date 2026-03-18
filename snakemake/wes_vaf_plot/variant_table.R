library(tidyverse)
library(gt)
library(glue)
library(here)

SOURCE_COL <- paste0("INFO_", snakemake@config$source_tag %||% "SOURCE")

format_for_gt <- function(tb, wanted_genes = NULL) {
  if (is.null(wanted_genes)) {
    wanted_genes <- tb$SYMBOL
  }
  tb |>
    filter(SYMBOL %in% wanted_genes) |>
    group_by(subject, SYMBOL, HGVSg) |>
    summarise(
      AD = first(AD),
      AF = mean(AF, na.rm = TRUE),
      DP = mean(INFO_DP, na.rm = TRUE),
      GT = first(GT),
      n_called_by = length(unique(!!as.symbol(SOURCE_COL))),
      Consequence = first(Consequence),
      Impact = first(IMPACT),
      HGVSp = str_remove(first(HGVSp), "ENS.*\\..*:p\\."),
      Consequence = first(Consequence),
    ) |>
    mutate(
      Consequence = str_replace_all(Consequence, "&", "<br>"),
      Consequence = str_replace_all(Consequence, "_", " "),
      Consequence = str_replace_all(Consequence, " variant$", ""),
      HGVSp = str_replace_all(HGVSp, "%3D", "=")
    ) |>
    ungroup()
}

make_table <- function(tb) {
  gt(
    tb,
    rowname_col = c("subject", "SYMBOL", "HGVSg"),
    row_group_as_column = TRUE,
  ) |>
    fmt_markdown(columns = "Consequence") |>
    fmt_number(columns = "AF") |>
    cols_label(
      "HGVSp" ~ "Amino acid change",
      "n_called_by" ~ "Number of callers",
      "AD" ~ "Alt. depth",
      "GT" ~ "Genotype",
      "AF" ~ "Alt. fraction",
      "DP" ~ "Depth"
    ) |>
    tab_stubhead(
      label = md(c("**Subject**", "**Gene**", "**HGVSp**"))
    ) |>
    cols_align("center", columns = c("GT", "n_called_by")) |>
    cols_label_with(fn = \(x) md(glue("**{x}**")))
}


variant_table <- function() {
  tb <- read_tsv(snakemake@input[[1]])
  rconfig <- snakemake@config[[snakemake@rule]]
  wanted_genes <- rconfig[["wanted_genes"]]
  format_for_gt(tb, wanted_genes) |>
    make_table() |>
    opt_table_font(size = 10) |>
    gtsave(snakemake@output[[1]])
}

variant_summary <- function() {
  tb <- read_tsv(snakemake@input[[1]]) |>
    mutate(HGVSg = str_remove(HGVSg, "^chr"))
  to_summarize <- snakemake@config$vars_to_summarize

  n_samples <- length(unique(tb$subject))

  kept <- lapply(names(to_summarize), \(symbol) {
    filtered <- tb |> filter(SYMBOL == symbol)
    spec <- to_summarize[[symbol]]
    lapply(names(spec), \(column) {
      filtered |> filter(!!as.symbol(column) %in% spec[[column]])
    }) |>
      bind_rows()
  }) |>
    bind_rows()

  prop_tb <- kept |>
    distinct(subject, HGVSg) |>
    group_by(HGVSg) |>
    count() |>
    mutate(prop = n / n_samples)

  fmtted <- kept |>
    distinct(subject, HGVSg, !!as.symbol(SOURCE_COL), .keep_all = TRUE) |>
    group_by(HGVSg) |>
    summarise(
      n_callers = map_dbl(list(unique(!!as.symbol(SOURCE_COL))), length),
      mean_af = mean(AF, na.rm = TRUE),
      GT = paste0(unique(discard(GT, \(x) x == ".")), collapse = ", "),
      SYMBOL = first(SYMBOL),
      HGVSp = first(HGVSp),
      cds_pos = first(CDS_position),
      Consequence = paste0(unique(Consequence), collapse = "&"),
    ) |>
    inner_join(prop_tb, by = join_by(HGVSg)) |>
    mutate(
      HGVSp = ifelse(is.na(HGVSp), NA, str_remove(HGVSp, ".*:")),
      HGVSg = ifelse(is.na(HGVSg), NA, str_remove(HGVSg, ".*:")),
    )

  tab <- fmtted |>
    gt(groupname_col = "SYMBOL", row_group_as_column = TRUE) |>
    fmt_number(columns = c("prop", "mean_af")) |>
    cols_label(
      n_callers = "N detecting callers",
      mean_af = "Mean AF",
      n = "N samples with variant",
      cds_pos = "CDS position",
      prop = "Proportion in cohort",
      GT = "Genotypes"
    ) |>
    cols_move_to_start(columns = c(HGVSg, HGVSp, cds_pos)) |>
    opt_stylize(style = 2) |>
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_row_groups(groups = everything())
    )
  gtsave(tab, snakemake@output[[1]])
}

if (exists("snakemake")) {
  globalenv()[[snakemake@rule]]()
} else {
  wanted_genes <- c("KRAS", "TP53")
  tb <- read_tsv(here(
    "2026-01-28_CRC",
    "vcfs",
    "combined.tsv"
  ))
  tab <- format_for_gt(tb, wanted_genes) |> make_table()
  tab |>
    opt_table_font(size = 10) |>
    gtsave("~/Downloads/variant_table.pdf")
}
