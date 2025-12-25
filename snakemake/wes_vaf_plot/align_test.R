suppressMessages({
  library(GenomicAlignments)
  library(glue)
  library(tidyverse)
  library(ggplot2)
  library(logger)
  library(ggmsa)
  library(here)
  library(patchwork)
})

reference <- readDNAStringSet(
  "/home/shannc/Bio_SDD/chula-stem/analyses/data_all/reference/genomes/GRCh38.p14_filtered.fasta"
)

chr <- 12
start <- 25245311
end <- 25245371
highlight <- 25245351
t_bam <- "/home/shannc/Bio_SDD/chula-stem/analyses/data_all/output/CCA/EXOME-PAIRED/13/tumor/4-13_tumor-recal.bam"
n_bam <- "/home/shannc/Bio_SDD/chula-stem/analyses/data_all/output/CCA/EXOME-PAIRED/13/normal/4-13_normal-recal.bam"

phred_ascii2pr <- function(ascii) {
  1 - 10^(-phred2ASCIIOffset(ascii) / 10)
}

make_pileup <- function(
  alignments,
  region,
  reference,
  prefix,
  with_prob = TRUE,
  target = NULL
) {
  if (length(region) != 1 || !"GRanges" %in% class(region)) {
    stop("`region` must be a GRanges object of length 1")
  }

  padding_char <- "+"
  piled <- GenomicAlignments::stackStringsFromGAlignments(
    alignments,
    region,
    Lpadding.letter = padding_char,
    Rpadding.letter = padding_char
  )
  to_remove <- as.character(XVector::subseq(piled, start = 1, width = 1)) ==
    padding_char
  piled <- piled[-to_remove]
  if (with_prob) {
    probs <- XVector::subseq(mcols(piled)$qual, start = target, width = 1) |>
      as.character() |>
      phred_ascii2pr() |>
      round(3)
    probs <- probs * 100 # Probability of correct base call at variant region
    align_names <- glue("{paste0(prefix, seq_along(piled))} ({probs})")
  } else {
    align_names <- paste0(prefix, seq_along(piled))
  }

  mcols(piled) <- NULL
  names(piled) <- align_names
  if (!is.null(reference)) {
    ref_subset <- Biostrings::DNAStringSet(XVector::subseq(
      reference[[seqnames(region)]],
      start = start(region),
      end = end(region)
    ))
    names(ref_subset) <- "REF"
    piled <- append(ref_subset, piled)
  }

  piled
}

region <- GenomicRanges::GRanges(
  seqnames = chr,
  ranges = IRanges::IRanges(start = start, end = end)
)

alignments <- lapply(c(t_bam, n_bam), \(bfile) {
  GenomicAlignments::readGAlignments(
    bfile,
    param = ScanBamParam(what = c("seq", "flag", "qual"), which = region)
  )
}) |>
  `names<-`(c("tumor", "normal"))


t_piled <- make_pileup(
  alignments$tumor,
  region,
  reference,
  "t",
  with_prob = TRUE,
  target = highlight - start + 1
)

n_piled <- make_pileup(
  alignments$normal,
  region,
  reference,
  "n",
  with_prob = TRUE,
  highlight - start + 1
)

make_plot <- function(pileup) {
  ggmsa(
    pileup,
    ref = "REF",
    consensus_views = TRUE,
    border = "white",
    font = "mono",
    position_highlight = highlight - start + 1,
    seq_name = TRUE,
    show.legend = FALSE
  ) +
    guides(fill = "none")
}


nplot <- make_plot(n_piled) + xlab("Normal")
tplot <- make_plot(t_piled) + xlab("Tumor") + scale_x_discrete(position = "top")

plot <- tplot +
  nplot +
  plot_annotation(
    title = "KRAS",
    theme = theme(
      plot.margin = unit(c(0, 0, 0, 0), "mm")
    )
  )
ggsave(here("test_msa.png"), plot, height = 10, width = 10)
