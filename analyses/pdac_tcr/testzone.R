library(Biostrings)
library(here)

testzone <- here("analyses", "output", "pdac_tcr", "testzone")
file <- here(testzone, "alignments_inter/10825037_TRA_1.fasta")

seqs <- readBStringSet(file)
subject <- seqs[names(seqs) == "sequence"]
patterns <- seqs[names(seqs) != "sequence"]

call <- list(
  pattern = patterns,
  subject = subject,
  type = "global-local",
  gapOpening = 50
)
align <- do.call(pwalign::pairwiseAlignment, call)
aligned_patterns <- pwalign::aligned(align)
names(aligned_patterns) <- names(patterns)
string_set <- c(subject, aligned_patterns)
writeXStringSet(string_set, here(testzone, "test_align.fasta"))

lapply()
