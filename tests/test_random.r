library(tidyverse)
library(GenomicRanges)
classify_cnv <- (
  "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/8-null-ClassifyCNV.tsv"
)
ref <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/aggregated_cnv.tsv"
mis <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/aggregated_msi.tsv"
out <- "/home/shannc/Bio_SDD/chula-stem/tests/msisensor/out.tsv"
dosage <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/Clingen-Dosage-Sensitivity-2024-12-16.csv"
clin <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/Clingen-Gene-Disease-Summary-2024-12-19.csv"
