library(tidyverse)
library(here)
library(glue)
source(here("analyses", "pdac", "main.R"))
vep_vcfs <- list.files(data_path, pattern = "7-P[0-9_]+-VEP_small.vcf.gz$", recursive = TRUE, full.names = TRUE)
bams <- list.files(data_path, pattern = "4-P[0-9_]+_tumor-recal.bam$", recursive = TRUE, full.names = TRUE)

wanted <- c(
  KRAS = "12:25205246-25250936",
  TP53 = "17:7661779-7687546",
  MUC5B = "11:1223066-1262172",
  KMT2C = "7:152134922-152436644",
  ARID1A = "1:26693236-26782104",
  SMAD4 = "18:51028528-51085045",
  GLI3 = "7:41960949-42264100",
  CDKN2A = "9:21967752-21995301"
)
target_file <- here("analyses", "pdac", "target_gene_regions.txt")

get_wanted_genes <- function(file) {
  args <- c("view", glue("--targets-file {target_file}"), file)
  pref <- utils$basename_no_ext(file) |> str_remove(".vcf")
  outfile <- glue("{pref}.vcf")
  filtered <- here(outdir, "vcfs", outfile)
  system2("bcftools", args, stdout = filtered)

  query_str <- "[%AD\t%AF\t%GT\t%PS]\t%INFO/DP\t%INFO/SOURCE"
  final_file <- here(outdir, "vcfs", glue("{pref}.tsv"))
  args2 <- c(glue("-i {filtered}"), glue("-o {final_file}"), glue("-q '{query_str}'"))
  system2(here("src", "bash", "get_vep_anno.bash"), args = args2)
}

get_wanted_bam <- function(file) {
  args <- c("view", "-b", "-h", glue("-L {target_file}"), file)
  pref <- utils$basename_no_ext(file)
  bam_out <- here(outdir, "bams", glue("{pref}.bam"))
  system2("samtools", args, stdout = bam_out)
}

if (path.expand("~") != "/home/shannc") {
  tmp <- lapply(vep_vcfs, get_wanted_genes)
  ## tmp <- lapply(bams, get_wanted_bam)
}

all_vcfs <- list.files(here(outdir, "vcfs"), pattern = "*.tsv", full.names = TRUE)
subjects <- utils$basename_no_ext(all_vcfs) |> str_extract("-(P.*)-VEP_small.*", group = 1)

combined_anno <- lapply(all_vcfs, \(x) {
  read_tsv(x) |>
    mutate(
      Protein_position = as.character(Protein_position),
      AF = as.double(AF),
      PS = as.double(PS),
      INFO_DP = as.double(INFO_DP)
    )
}) |>
  `names<-`(subjects) |>
  list_rbind(names_to = "subject")
combined_anno |> write_csv(here(outdir, "select_genes_VEP.csv"))

for (n in names(wanted)) {
  combined_anno |>
    filter(SYMBOL == n) |>
    group_by(subject, HGVSg) |>
    select(subject, HGVSg, HGVSp, GT, AD, AF, INFO_DP, INFO_SOURCE, Amino_acids, Existing_variation) |>
    arrange(subject, HGVSg) |>
    distinct() |>
    write_tsv(here(outdir, "formatted_tsv", glue("{n}.tsv")))
}


combined_anno |> distinct() |>
