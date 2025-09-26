library(tidyverse)
library(here)

# [2025-08-26 Tue] Aggregating results for select variants
vep_file <- here(
  "analyses",
  "output",
  "pdac_phcase",
  "8-PHCase_6-VEP_small.tsv"
)
if (!file.exists(vep_file)) {
  file.copy(
    here(
      "output",
      "PHCase",
      "PHCase_6",
      "annotations",
      "8-PHCase_6-VEP_small.tsv"
    ),
    vep_file
  )
} else {
  vep <- read_tsv(vep_file)
}

wanted_symbols <- c("ARID1A", "BRAF", "TET2", "MAP2K4", "TP53")

# TODO: need a way to choose a single transcript only, to make reporting easier
#

# Agreement with liquid biopsy variants
# ARID1A: Yes
# BRAF: Yes, ENSP00000419060.2:p.Asn486_Pro490del
# MAP2K4: Yes
# TP53: Yes
# TET2: not found

filtered <- vep |>
  filter(SYMBOL %in% wanted_symbols)
