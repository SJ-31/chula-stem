library(cnp.mops)
library(here)

normal <- here("analyses", "data_all", "output", "WES_PON", "HCC_19", "normal", "4-HCC_19_normal-recal.bam")
tumor <- here("analyses", "data_all", "output", "HCC", "Exome", "P17", "tumor", "4-P17_tumor-recal.bam")

n_ranges <- getReadCountsFromBam(normal)
t_ranges <- getReadCountsFromBam(tumor)

# This doesn't seem compatible with the exome variant...
result <- referencecn.mops(t_ranges, n_ranges)

## Use referencn.mops()
