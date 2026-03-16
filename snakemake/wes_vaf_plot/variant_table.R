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
