library(tidyverse)
library(here)
library(glue)

rnaseq_dir <- here("analyses", "data_all", "output", "PDAC-RNASEQ")
phcase_13 <- here(rnaseq_dir, "PHcase_13", "tumor")
# phcase13 has avg depth of 4.30
phcase_14 <- here(rnaseq_dir, "PHcase_14", "tumor")
# phcase14 has avg depth of 5.13

get_read_counts <- function(sample, dir) {
  read_tsv(
    here(dir, glue("2-{sample}_tumor-STAR_Counts.tsv")),
    col_names = c("id", "count")
  )
}

counts <- get_read_counts("PHcase_13", phcase_13) |>
  inner_join(
    get_read_counts("PHcase_14", phcase_14),
    by = "id",
    suffix = c("_13", "_14")
  )
avg_diff <- mean(counts$count_13 - counts$count_14)

coverage <- inner_join(
  read_tsv(here(phcase_13, "coverage.txt")),
  read_tsv(here(phcase_14, "coverage.txt")),
  by = "#rname",
  suffix = c("_13", "_14")
)

compare_cov <- coverage |>
  select("#rname", contains("coverage"), contains("meandepth")) |>
  filter(grepl("^([0-9XY]+|MT)", `#rname`))
# coverage of PHcase_14 is nearly twice that of PHcase_13
