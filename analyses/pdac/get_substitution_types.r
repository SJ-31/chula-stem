library(MutationalPatterns)
library(tidyverse)
library(here)
here::i_am("./analyses/pdac/vaf_plot.r")

sbs_merged_file <- here("analyses", "output", "pdac_sbs_mp.tsv")
data_path <- here("analyses", "data_all", "output", "PDAC")

if (!file.exists(sbs_merged_file)) {
  vcfs <- list.files(data_path, pattern = "7-P[0-9_]+-Small_high_conf.vcf.gz", recursive = TRUE, full.names = TRUE)
  ref <- "BSgenome.Hsapiens.UCSC.hg38"
  library(ref, character.only = TRUE)
  sample_names <- basename(vcfs) |> str_extract("P[0-9_]+")
  granges <- MutationalPatterns::read_vcfs_as_granges(vcfs, sample_names = sample_names, ref)
  sbs <- lapply(sample_names, \(x) {
    gr <- granges[[x]]
    muts <- mut_type(gr) |> discard(\(x) str_length(x) != 3)
    tibble(type = muts) |>
      group_by(type) |>
      summarize(count = n()) |>
      mutate(sample = x)
  }) |> bind_rows()
  write_tsv(sbs, sbs_merged_file)
}

# <2025-01-23 Thu> the results are pretty much identical to the substitution values
# obtained with bcftools stats, so just use that
