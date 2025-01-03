library(tidyverse)
library(GenomicRanges)
library(glue)
library(grid)
library(ggplot2)
library(here)
src <- Sys.getenv("stem_r_src")
source("/home/shannc/Bio_SDD/chula-stem/src/R/plotting.R")
## source(glue("{src}/plotting.R"))

cns_file <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.call.cns"
# Contains copy number calls
cnr_file <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.cnr"
# Contains the log2 segmentation data
cnn_file <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.targetcoverage.cnn"
targets <- ""

facets_rds <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/5-patient_10-Facets/5-patient_10-Facets_hisens.rds"

facets <- readRDS(facets_rds)
# Off-target bins typically going to be much wider than target bins

regions_bed <- "/home/shannc/Bio_SDD/chula-stem/tests/data/Unzipped_regions.bed"
