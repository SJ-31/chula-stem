library(ggmsa)
library(ggplot2)
library(here)
library(Biostrings)
library(patchwork)

workdir <- here("analyses", "marmoset")
# MANE select transcript for human is NM_000377.3

prot <- here(workdir, "wasp_prot_aligned.fasta")
nuc <- here(workdir, "wasp_nuc_aligned.fasta")

prot_set <- readAAMultipleAlignment(prot)

interval <- 100

plot_msa_multiple <- function(file, width, seqtype = "prot", ...) {
  library(patchwork)
  if (seqtype == "prot") {
    obj <- Biostrings::readAAMultipleAlignment(file)
  } else {
    obj <- Biostrings::readDNAMultipleAlignment(file)
  }
  len <- nchar(obj)
  sections <- seq(1, len, by = width)
  last <- tail(sections, n = 1)
  if ((len - last) < (width / 2)) {
    sections[length(sections)] <- len
  } else if (tail(sections, n = 1) != len) {
    sections <- c(sections, len)
  }
  last <- tail(sections, n = 1)
  plots <- lapply(seq_along(sections), \(i) {
    if (sections[i] == last) {
      NULL
    } else {
      ggmsa::ggmsa(file, sections[i], sections[i + 1], ...)
    }
  })
  purrr::reduce(plots, \(x, y) x / y)
}

prot_plot <- plot_msa_multiple(prot, 100, seq_name = TRUE)
nuc_plot <- plot_msa_multiple(nuc, 200, seq_name = TRUE, seqtype = "nuc")

ggsave(
  here(workdir, "protein_aligned.png"),
  prot_plot,
  width = 30,
  dpi = 500
)

ggsave(
  here(workdir, "nuc_aligned.png"),
  nuc_plot,
  width = 30,
  height = 10,
  dpi = 500
)
