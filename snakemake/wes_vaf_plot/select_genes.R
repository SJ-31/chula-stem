library(tidyverse)
library(here)
library(glue)
source(here("src", "R", "utils.R"))
rlang::global_entrace(enable = TRUE)

data_mapping <- snakemake@params$data_mapping
cases <- snakemake@params$cases

SOURCE_TAG <- snakemake@config$source_tag
if (is.null(SOURCE_TAG)) {
  SOURCE_TAG <- "TOOL_SOURCE"
}
OUTDIR <- snakemake@params$outdir
TARGET_FILE <- glue("{snakemake@scriptdir}/target_genes.tsv")


vcfs <- sapply(
  cases,
  \(x) {
    glue("{data_mapping[[x]]}/{x}/annotations/7-{x}-VEP_small.vcf.gz")
  },
  simplify = FALSE,
  USE.NAMES = TRUE
)


gene_tb <- unlist(snakemake@params$genes, use.names = FALSE) |>
  lapply(\(x) {
    splits <- str_split_1(x, ":")
    chrom <- splits[1]
    start_end <- str_split_1(splits[2], "-")
    tibble(chrom = chrom, start = start_end[1], end = start_end[2])
  }) |>
  bind_rows()
write_tsv(gene_tb, TARGET_FILE, col_names = FALSE)

get_wanted_genes <- function(name, file, is_paired = TRUE, vep = TRUE) {
  args <- c(
    "view",
    glue("--targets-file {TARGET_FILE}")
  )
  if (is_paired) {
    samples <- system2("bcftools", args = c("query", "-l", file), stdout = TRUE)
    tumor_sample <- keep(samples, \(x) str_detect(x, "tumor"))[1]
    args <- c(args, glue("-s {tumor_sample}"))
  }
  args <- c(args, file)
  tmp <- glue("{OUTDIR}/{name}.vcf")

  system2("bcftools", args, stdout = tmp)

  outfile <- glue("{OUTDIR}/{name}.tsv")
  query <- glue(
    "[%AD\t%AF\t%GT\t%PS]\t%INFO/DP\t%INFO/{SOURCE_TAG}\t%CHROM\t%POS"
  )
  args2 <- c(
    glue("-i {tmp}"),
    glue("-o {outfile}"),
    glue("-q '{query}'")
  )
  if (vep) {
    system2(here("src", "bash", "get_vep_anno.bash"), args = args2)
  }
  read_tsv(outfile)
}

vep_vcfs <- lapply(names(vcfs), \(x) {
  get_wanted_genes(
    name = x,
    file = vcfs[[x]],
    vep = TRUE,
    is_paired = snakemake@config$variant_calling$paired
  )
}) |>
  `names<-`(names(vcfs))


to_double <- c("AF", "PS", "INFO_DP")

combined_anno <- lapply(
  vep_vcfs,
  \(x) {
    cns <- colnames(x)
    to_character <- cns[!cns %in% to_double]
    mutate(x, across(all_of(to_double), as.double)) |>
      mutate(across(all_of(to_character), as.character))
  }
) |>
  list_rbind(names_to = "subject")

write_tsv(combined_anno, snakemake@output$combined)

unlink(TARGET_FILE)
