library(tidyverse)
library(glue)
library(ggplot2)
library(here)

msi <- read_tsv(here(
  "analyses",
  "output",
  "crc_exome",
  "6-CEN2-Msisensor_summaries_paired.tsv"
))
sbs <- read_tsv(here(
  "analyses",
  "output",
  "crc_exome",
  "8-CEN2-SigProfilerAssignment_activities_paired.tsv"
))


## ggplot(msi, aes(x = Total_Number_of_Sites, y = `%`)) +
##   geom_point()

sbs_filtered <- sbs |>
  select(where(\(x) {
    if (is.character(x)) {
      TRUE
    } else {
      sum(x) > 100
    }
  })) # [2025-08-06 Wed] Got SBS1 and SBS5, which aren't particularly interesting
