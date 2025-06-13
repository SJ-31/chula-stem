library(here)
library(GenomicAlignments)
library(tidyverse)
library(ggplot2)
library(ggmsa)

## * Setup
ref_genome <- here("analyses", "data", "GRCh38.p14_filtered.fasta")
outdir <- here("analyses", "output", "pdac")
bamdir <- here(outdir, "bams")

testfile <- here(bamdir, "4-P1_tumor-recal.bam")
read <- GenomicAlignments::readGAlignments(
  testfile,
  param = ScanBamParam(what = c("seq", "flag", "qual"))
)
seqlevels(read) <- paste0("chr", seqlevels(read))
read <- read[seqnames(read) == "chr12", ]
ranges <- granges(read)

genome <- readRDS(here("analyses", "data", "GRCh38.p14_filtered.rds"))
seqlevels(genome) <- paste0("chr", seqlevels(genome))


## * Plot
# %%

#' Plot alignments in a BAM file with ggmsa
#'
#' @description
#' @param galignment A GAlignments object from the "GenomicAlignments" package
#'    Must have "seq" in its mcols, which you can get by passing
#'   param = ScanBamParam(what = c("seq")) to the readGAlignments call
#' @param start Start index of region to plot
#' @param end End index of region to plot
#' @param seq Name of sequence
#' @param reference Reference genome used to make the alignments
plot_alignment_as_msa <- function(
  galignment,
  start,
  end,
  seq,
  reference = NULL,
  prefix = "a",
  ...
) {
  region <- GRanges(
    seqnames = seq,
    ranges = IRanges(start = start, end = end)
  )
  piled <- stackStringsFromGAlignments(galignment, region = region)
  mcols(piled) <- NULL
  names(piled) <- paste0(prefix, seq_along(piled))
  if (!is.null(reference)) {
    stopifnot(
      "Reference must be a Biostrings object!" = inherits(
        reference,
        "XStringSet"
      ),
      "Given sequence not in reference!" = seq %in% names(reference)
    )
    ref_seq <- DNAStringSet(subseq(reference[[seq]], start = start, end = end))
    names(ref_seq) <- "reference"
    piled <- append(ref_seq, piled)
    ref_str <- "reference"
  } else {
    ref_str <- NULL
  }
  print(piled)
  ggmsa::ggmsa(
    piled,
    ref = ref_str,
    consensus_views = !is.null(ref_str),
    ...
  )
}

chrom <- "chr12"
from <- 25245345
to <- from + 10

msa <- plot_alignment_as_msa(
  read,
  start = from,
  end = to,
  seq = chrom,
  reference = genome,
  font = "mono",
  color = "Chemistry_NT"
) +
  labs(title = "")

msa +
  scale_x_continuous(breaks = seq(from, to, by = 3)) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12, angle = -90)
  ) %>%
    ggsave(plot = ., filename = here(outdir, "test.png"))
