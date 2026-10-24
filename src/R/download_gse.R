library(tidyverse)
library(here)
library(glue)
library(GEOquery)

if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(
    parser,
    c("-i", "--input"),
    type = "character",
    help = "Input filename or accession"
  )

  args <- parse_args(parser)
  is_file <- if (file.exists(args$input)) TRUE else FALSE
  if (is_file) {
    gse <- getGEO(GEO = NULL, filename = args$input)
  } else {
    gse <- getGEO(args$input)
  }
  if (length(gse) == 1) {
    pdata <- if ("ExpressionSet" %in% class(gse)) {
      pData(gse)
    } else {
      pData(gse[[1]])
    }
    as_tibble(pdata) |>
      write_tsv(here("meta", glue("{args$input}.tsv")))
  } else {
    for (n in names(gse)) {
      if (str_detect(n, "_series_matrix")) {
        name <- str_remove(n, "_series_matrix.*")
        pData(gse[[n]]) |>
          as_tibble() |>
          write_tsv(here("meta", glue("{name}.tsv")))
      }
    }
  }
}
