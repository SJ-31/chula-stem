library(tidyverse)
library(GenomicRanges)
library(glue)
library(grid)
library(ggplot2)
library(here)
here::i_am("tests/test_random.r")
src <- Sys.getenv("stem_r_src")
source(here("src", "R", "utils.R"))
## source(glue("{src}/plotting.R"))

cns_file <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.call.cns"
# Contains copy number calls
cnr_file <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.cnr"
# Contains the log2 segmentation data
vep_small_file <- here("tests", "vep", "7-patient_10-VEP_small_1.tsv")
vep_sv_file <- here("tests", "vep", "sv.tsv")
cnn_file <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/4-patient_10_cancer-recal.targetcoverage.cnn"
targets <- ""
classify_cnv <- here("tests", "classify_cnv", "4-classify_cnv-CR.tsv")

facets_rds <- "/home/shannc/Bio_SDD/chula-stem/tests/classify_cnv/5-patient_10-Facets/5-patient_10-Facets_hisens.rds"

facets <- readRDS(facets_rds)
# Off-target bins typically going to be much wider than target bins

regions_bed <- "/home/shannc/Bio_SDD/chula-stem/tests/data/Unzipped_regions.bed"


## * Setup for plotting cn data <2025-01-10 Fri>
library(Gviz)
library(ensembldb)
library(AnnotationHub)

source(here("src", "R", "plotting.R"))
setAnnotationHubOption("CACHE", here(".cache", "AnnotationHub"))
ah <- AnnotationHub()
obj <- "AH116860" # Ensembl 112 EnsDb for Homo sapiens
ensdb <- ah[[obj]]
seqlevelsStyle(ensdb) <- "UCSC"

ratios <- read_tsv(cnr_file) |> format_chr_data(divide_by = 1)
cn_calls <- read_tsv(cns_file) |> format_chr_data(divide_by = 1)
ctb <- read_tsv(classify_cnv) |> rename_with(\(x) str_replace_all(x, "[ -]", "_"))

chromosome_from_loc <- function(tb, loc_col = "Loc") {
  separated <- separate_wider_delim(tb, all_of(loc_col),
    delim = ":", names = c("chromosome", "start")
  )
  if (!str_detect(separated$chromosome[1], "chr")) {
    separated$chromosome <- paste0("chr", separated$chromosome)
  }
  separated |>
    mutate(
      start = as.numeric(start),
      end = start + (nchar(Alt) - nchar(Ref)),
      end = ifelse(end - start <= 0, start + 1, end)
    ) |>
    relocate(chromosome, start, end, .before = everything())
}

get_cn_iranges <- function(cn_tb, chr) {
  filtered <- dplyr::filter(cn_tb, chromosome == chr)
  with(filtered, IRanges(start = start, end = end, names = cn))
}

vep_svs <- read_tsv(vep_sv_file, col_types = "c") |>
  chromosome_from_loc() |>
  mutate(VAF = as.character(VAF))
vep_small <- read_tsv(vep_small_file, col_types = "c") |> chromosome_from_loc()
all_vep <- bind_rows(vep_svs, vep_small)

joined <- ratios |>
  inner_join(cn_calls,
    by = join_by(chromosome, between(mid, y$start, y$end)),
    suffix = c("", "_right")
  )

## ** Plots

## *** Get granges objects
copy_ratios_granges <- makeGRangesFromDataFrame(ratios)
cn_granges <- makeGRangesFromDataFrame(cn_calls)
var_granges <- makeGRangesFromDataFrame(all_vep, strand.field = "*")
mcols(var_granges)$count <- 1
mcols(copy_ratios_granges)$log2 <- ratios$log2
mcols(cn_granges)$cn <- cn_calls$cn

genes <- list()
current_chr <- "chr9"
types <- c("ds", "all")
cols <- c("Known_or_predicted_dosage_sensitive_genes", "All_protein_coding_genes")
for (i in seq(2)) {
  genes[[types[i]]] <- ctb |>
    dplyr::filter((Chromosome == current_chr) &
      (!is.na(all_of(cols[i])))) |>
    pluck(cols[i]) |>
    flatten_by(separator = ",", collapse = FALSE) |>
    unlist() |>
    map_chr(trimws) |>
    unique()
}

all_ctb <- ctb$All_protein_coding_genes |>
  flatten_by(",", FALSE) |>
  unlist() |>
  map_chr(trimws) |>
  unique()

# Filter genes on the current chromosome and those flagged as dosage-sensitive
## genes_filtered <- genes(ensdb, filter = ~ gene_name %in% genes$all & seq_name == current_chr)

## all_vep |>
##   dplyr::filter(SYMBOL %in% genes_filtered$symbol) |>
##   pluck("SYMBOL")

## ** Filter on current gene
# This could be replaced with a range as well
current_gene <- genes_filtered[81, ]

# Get features overlapping with the gene
gene_data <- getGeneRegionTrackForGviz(ensdb, filter = GRangesFilter(current_gene))
var_filtered <- subsetByOverlaps(var_granges, current_gene)
cn_filtered <- pintersect(cn_granges, current_gene, drop.nohit.ranges = TRUE)
# TODO: what is the difference between pintersect and intersect???
# pintersect performs an

copy_ratios_filtered <- subsetByOverlaps(copy_ratios_granges, current_gene)
copy_ratios_filtered <- pintersect(copy_ratios_granges, current_gene, drop.nohit.ranges = TRUE)
# <2025-01-10 Fri> pintersect is better for this, but why???

# Make tracks
gtrack <- GenomeAxisTrack(cn_filtered)
copy_data_track <- DataTrack(copy_ratios_filtered, name = "Copy Ratio")
var_track <- DataTrack(var_filtered, type = "hist", name = "Variant count")
gene_track <- GeneRegionTrack(gene_data, name = "Features")

tracks <- list(gtrack, gene_track, copy_data_track, var_track)
# <2025-01-10 Fri> Plan is just to pass it the cnv regions exactly,
# and extend using the plot track params
# Want three tracks: ideaogram, gene models, and cn data

plotTracks(tracks,
  transcriptAnnotation = "transcript", collapseTranscripts = TRUE,
  shape = "arrow"
)
