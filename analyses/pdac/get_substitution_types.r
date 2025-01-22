library(MutationalPatterns)
here::i_am("./analyses/pdac/vaf_plot.r")
source(here("src", "R", "utils.R"))

sbs_merged_file <- here("analyses", "output", "pdac_sbs_all.tsv")
data_path <- here("analyses", "data_all", "output", "PDAC")

vcfs <- list.files(data_path, pattern = "7-P[0-9_]+-Small_high_conf.vcf.gz", recursive = TRUE, full.names = TRUE)

sample <- vcfs[1]
gr <- MutationalPatterns::read_vcfs_as_granges(sample, sample_names = "P1", )
