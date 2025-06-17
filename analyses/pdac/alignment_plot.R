library(here)
library(GenomicAlignments)
library(glue)
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

genome <- readRDS(here("analyses", "data", "GRCh38.p14_filtered.rds"))
seqlevels(genome) <- paste0("chr", seqlevels(genome))


## * Plot
# %%

#' Plot alignments in a BAM file with ggmsa
#'
#' TODO: The position highlight could be better
#' NOTE: the pileup formatting is good, but the plot is kinda meh, only useful
#'    for a single snp
#' you could improve this by replacing ggmsa with a custom plotting engine
#' @description
#' @param galignment A GAlignments object from the "GenomicAlignments" package
#'    Must have "seq" in its mcols, which you can get by passing
#'   param = ScanBamParam(what = c("seq")) to the readGAlignments call
#' @param start Start index of region to plot
#' @param end End index of region to plot
#' @param seq Name of sequence
#' @param reference Reference genome used to make the alignments, as a Biostrings object
#'    e.g. DNAString set
#' @param position_highlight position to highlight by ggmsa
#' @param ... Additional arguments passed to ggmsa
plot_alignment_as_msa <- function(
  galignment,
  start,
  end,
  seq,
  reference = NULL,
  prefix = "a",
  position_highlight = NULL,
  ...
) {
  region <- GenomicRanges::GRanges(
    seqnames = seq,
    ranges = IRanges::IRanges(start = start, end = end)
  )
  piled <- GenomicAlignments::stackStringsFromGAlignments(
    galignment,
    region = region
  )
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
    ref_seq <- Biostrings::DNAStringSet(XVector::subseq(
      reference[[seq]],
      start = start,
      end = end
    ))
    names(ref_seq) <- "reference"
    piled <- append(ref_seq, piled)
    ref_str <- "reference"
  } else {
    ref_str <- NULL
  }
  if (!is.null(position_highlight)) {
    position_highlight <- position_highlight - from + 1
  }
  ggmsa::ggmsa(
    piled,
    ref = ref_str,
    consensus_views = !is.null(ref_str),
    position_highlight = position_highlight,
    ...
  )
}


chrom <- "chr12"
from <- 25245345
to <- from + 10
highlight <- c(25245349, 25245350, 25245351)


plot_helper <- function(bam_file, extract = "4-(P[0-9_]+_tumor)-.*") {
  subject <- basename(bam_file) |>
    str_extract(extract, group = 1)
  bam <- GenomicAlignments::readGAlignments(
    bam_file,
    param = ScanBamParam(what = c("seq", "flag", "qual"))
  )
  bam <- bam[seqnames(bam) == "chr12", ]
  msa <- plot_alignment_as_msa(
    bam,
    start = from,
    end = to,
    seq = chrom,
    reference = genome,
    position_highlight = highlight,
    font = "mono",
    color = "Chemistry_NT",
  ) +
    labs(
      title = "KRAS Gly12 Alignments",
      subtitle = glue("Subject: {subject}")
    ) +
    guides(fill = "none")
  ggsave(
    plot = msa,
    filename = here(outdir, "alignment_png", glue("{subject}.png"))
  )
}

bams <- list.files(bamdir, pattern = "4-*.bam$", full.names = TRUE)
bams_old <- list.files(bamdir, pattern = "^N.*.bam$", full.names = TRUE)
tmp <- lapply(bams, plot_helper)
tmp <- lapply(bams_old, \(x) plot_helper(x, extract = "(N[0-9_BC]+).bam"))
