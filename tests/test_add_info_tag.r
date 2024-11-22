library("vcfR")
library("glue")

vcf_info_add_tag <- function(name, description, number, type, default, input, output) {
  vcf <- read.vcfR(input, verbose = FALSE)
  info <- vcf@fix[, "INFO"]
  vcf@fix[, "INFO"] <- paste0(info, glue(";SOURCE={default}"))
  m <- glue("##INFO=<ID={name},Number={number},Type={type},Description=\"{description}\">")
  vcf@meta <- append(vcf@meta, m)
  write.vcf(vcf, output)
}

## out <- "/home/shannc/Bio_SDD/chula-stem/tests/added.vcf.gz"
## vcf_info_add_tag("SOURCE", "blahblah", ".", "String", "mutect2", sample, out)
