library(GenomicAlignments)
library(glue)
library(tidyverse)
library(ggplot2)
library(ggmsa)


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

plot_helper <- function(
  subject,
  bam,
  start,
  end,
  chrom,
  ref_genome,
  title,
  outdir,
  highlight = NULL,
  extract = "4-(P[0-9_]+_tumor)-.*",
  suffix = ""
) {
  bam_current <- bam[seqnames(bam) == chrom, ]
  msa <- plot_alignment_as_msa(
    bam_current,
    start = start,
    end = end,
    seq = chrom,
    reference = ref_genome,
    position_highlight = highlight,
    font = "mono",
    color = "Chemistry_NT",
  ) +
    labs(
      title = title,
      subtitle = glue("Subject: {subject}")
    ) +
    guides(fill = "none")
  filename <- glue("{subject}_{chrom}:{start}-{end}:{suffix}")
  ggsave(
    plot = msa,
    filename = glue("{outdir}/{filename}.png")
  )
}

cfg <- snakemake@config$alignment_plot
ref_genome <- readDNAStringSet(cfg$reference_genome)
if (cfg$ref_add_chr) {
  seqlevels(ref_genome) <- paste0("chr", seqlevels(ref_genome))
}
outdir <- snakemake@params$outdir
data_path <- snakemake@params$data_path

samples <- snakemake@input
sample_tb <- tibble(vcf = samples) |>
  mutate(sample = str_remove(basename(vcf), ".tsv")) |>
  mutate(
    bam = map_chr(sample, \(s) {
      glue("{data_path}/{s}/tumor/4-{s}_tumor-recal.bam")
    })
  )

t0 <- apply(sample_tb, 1, \(tb_row) {
  bam_file <- tb_row["bam"]
  sample <- tb_row["sample"]
  vcf <- read_tsv(tb_row["vcf"])
  bam <- GenomicAlignments::readGAlignments(
    bam_file,
    param = ScanBamParam(what = c("seq", "flag", "qual"))
  )
  for (symbol in cfg$target_symbols) {
    filtered <- vcf |>
      filter(SYMBOL == symbol) |>
      select(CHROM, POS, HGVSg, HGVSp) |>
      distinct()
    t1 <- apply(filtered, 1, \(row) {
      chrom <- row["CHROM"]
      if (cfg$sample_add_chr) chrom <- glue("chr{chrom}")
      loc <- row["POS"]
      hgvsg <- row["HGVSg"]
      hgvsp <- row["HGVSp"]
      title <- glue("{symbol}  HGVSg: {hgvsg}")
      if (!is.null(hgvsp)) title <- glue("{title}, HGVSp: {hgvsp}")
      plot_helper(
        bam = bam,
        start = max(loc - cfg$window, 1),
        outdir = glue("{outdir}/{sample}"),
        end = loc + cfg$window,
        chrom = chrom,
        ref_genome = ref_genome,
        title = title,
        highlight = loc,
        suffix = symbol
      )
    })
  }
})
