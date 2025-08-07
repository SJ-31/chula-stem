library(dndscv)
library(readr)

tb <- read.table(snakemake@input[[1]])

get <- function(list, val, default) {
  v <- list[[val]]
  if (is.null(v)) {
    default
  } else {
    v
  }
}

config <- snakemake@config
output <- snakemake@output

if (!is.null(config$gene_list)) {
  gene_list <- read_lines(config$gene_list)
} else {
  gene_list <- NULL
}

if (!is.null(config$kc)) {
  kc <- read_lines(config$kc)
} else {
  kc <- "cgc81"
}

result <- dndscv(
  tb,
  refdb = get(config, "refdb", "hg38"),
  sm = get(config, "sm", "192r_3w"),
  kc = kc,
  max_muts_per_gene_per_sample = get(config, "max_muts_per_gene_per_sample", 3),
  max_coding_muts_per_sample = get(config, "max_coding_muts_per_sample", 3000),
  use_indel_sites = get(config, "use_indel_sites", TRUE),
  min_indels = get(config, "min_indels", 5),
  maxcovs = get(config, "maxcovs", 20),
  constrain_wnon_wspl = get(config, "constrain_wnon_wspl", TRUE),
  numcode = get(config, "numcode", 1),
  mingenecovs = get(config, "mingenecovs", 500),
  gene_list = gene_list
)

tsvs <- c(
  "globaldnds",
  "sel_cv",
  "sel_loc",
  "annotmuts",
  "genemuts",
  "geneindels",
  "mle_submodel",
  "exclmuts",
  "exclsamples",
  "wrongmuts"
)

rds <- c(
  "nbreg",
  "nbregind",
  "poissmodel"
)

for (tsv in tsvs) {
  write_tsv(result[[tsv]], output[[tsv]])
}

for (obj in rds) {
  saveRDS(result[[obj]], output[[obj]])
}

write_lines(result$exclsamples, output$exclsamples)
