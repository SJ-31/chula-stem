suppressMessages({
  library(GenomicAlignments)
  library(glue)
  library(tidyverse)
  library(ggplot2)
  library(logger)
  library(ggmsa)
  library(patchwork)

  log_threshold(TRACE)
  file <- appender_file(snakemake@log[[1]])
  log_appender(file)
})

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

  piled <- GenomicAlignments::stackStringsFromGAlignments(alignments, region)
  if (with_prob) {
    probs <- XVector::subseq(mcols(piled)$qual, start = target, end = target) |>
      as.character() |>
      phred_ascii2pr() |>
      round(3)
    probs <- probs * 100
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

#' Plot alignments in a BAM file with ggmsa
#'
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
#' @param reference Reference genome used to make the alignments,
#' as a Biostrings object    e.g. DNAString set
#' @param position_highlight position to highlight by ggmsa
#' @param ... Additional arguments passed to ggmsa
plot_alignment_as_msa <- function(
  bam_file,
  start,
  end,
  seq,
  title,
  subtitle,
  reference = NULL,
  normal_bam_file = "",
  position_highlight = NULL,
  ...
) {
  region <- GenomicRanges::GRanges(
    seqnames = seq,
    ranges = IRanges::IRanges(start = start, end = end)
  )
  alignments <- lapply(c(bam_file, normal_bam_file), \(bfile) {
    if (!file.exists(bfile)) {
      NULL
    } else {
      GenomicAlignments::readGAlignments(
        bfile,
        param = ScanBamParam(what = c("seq", "flag", "qual"), which = region)
      )
    }
  }) |>
    `names<-`(c("tumor", "normal"))

  tumor_piled <- make_pileup(alignments$tumor, region, reference, "t")
  if (!is.null(alignments$normal)) {
    normal_piled <- make_pileup(alignments$normal, region, reference, "n")
  } else {
    normal_piled <- NULL
  }
  if (!is.null(reference)) {
    stopifnot(
      "Reference must be a Biostrings object!" = inherits(
        reference,
        "XStringSet"
      ),
      "Given sequence not in reference!" = seq %in% names(reference)
    )
    ref_str <- "REF"
  } else {
    ref_str <- NULL
  }
  if (!is.null(position_highlight)) {
    position_highlight <- position_highlight - start + 1
  }

  make_plot <- function(pileup) {
    ggmsa(
      pileup,
      ref = ref_str,
      consensus_views = !is.null(ref_str),
      border = "white",
      position_highlight = position_highlight,
      show.legend = FALSE,
      ...
    ) +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank()
      ) +
      guides(fill = "none")
  }

  tplot <- make_plot(tumor_piled)
  if (is.null(normal_piled)) {
    tplot + labs(title = title, subtitle = subtitle)
  } else {
    tplot <- tplot + xlab("Tumor") + scale_x_discrete(position = "top")
    nplot <- make_plot(normal_piled) +
      xlab("Normal") +
      scale_x_discrete(position = "top")
    tplot +
      nplot +
      plot_annotation(
        title = title,
        subtitle = subtitle
      )
  }
}

plot_helper <- function(
  subject,
  bam_file,
  normal_bam_file,
  start,
  end,
  chrom,
  ref_genome,
  title,
  outdir,
  highlight = NULL,
  extract = "4-(P[0-9_]+_tumor)-.*",
  suffix = "",
  height = 12,
  width = 12
) {
  log_debug("Plotting for {start}-{end}")
  msa <- plot_alignment_as_msa(
    bam_file = bam_file,
    normal_bam_file = normal_bam_file,
    start = start,
    end = end,
    seq = chrom,
    reference = ref_genome,
    title = title,
    subtitle = glue("Subject: {subject}"),
    position_highlight = highlight,
    font = "mono",
    color = "Chemistry_NT",
  )
  filename <- glue("{subject}_{chrom}:{start}-{end}:{suffix}")
  ggsave(
    plot = msa,
    filename = glue("{outdir}/{filename}.png"),
    height = height,
    width = width
  )
  knitr::plot_crop(glue("{outdir}/{filename}.png"))
}

cfg <- snakemake@config$alignment_plot
ref_genome <- readDNAStringSet(cfg$reference_genome)
window <- as.numeric(cfg$window)
log_debug("window {window}")
outdir <- snakemake@params$outdir
data_path <- snakemake@params$data_path

samples <- snakemake@input$vcfs
sample_tb <- tibble(vcf = samples) |>
  mutate(sample = str_remove(basename(vcf), ".tsv")) |>
  mutate(
    bam = glue("{data_path}/{sample}/tumor/4-{sample}_tumor-recal.bam"),
    normal_bam = glue(
      "{data_path}/{sample}/normal/4-{sample}_normal-recal.bam"
    ),
  )

t0 <- apply(sample_tb, 1, \(tb_row) {
  bam_file <- tb_row["bam"]
  normal_bam_file <- tb_row["normal_bam"]
  sample <- tb_row["sample"]
  vcf <- read_tsv(tb_row["vcf"])
  dir.create(glue("{outdir}/{sample}"), showWarnings = FALSE)
  for (symbol in cfg$target_symbols) {
    filtered <- vcf |>
      filter(SYMBOL == symbol & !is.na(POS)) |>
      select(CHROM, POS, HGVSg, HGVSp) |>
      distinct()
    if (nrow(filtered) == 0) {
      next
    }
    t1 <- apply(filtered, 1, \(row) {
      log_debug("row: {row}")
      chrom <- row["CHROM"]
      loc <- as.numeric(row["POS"])
      hgvsg <- row["HGVSg"]
      hgvsp <- row["HGVSp"]
      title <- glue("{symbol}  HGVSg: {hgvsg}")
      log_debug("start: {max(loc - cfg$window, 1)}")
      log_debug("end: {loc + cfg$window}")
      if (!is.null(hgvsp)) title <- glue("{title}, HGVSp: {hgvsp}")
      plot_helper(
        subject = sample,
        bam_file = bam_file,
        normal_bam_file = normal_bam_file,
        start = max(loc - window, 1),
        end = loc + window,
        outdir = glue("{outdir}/{sample}"),
        chrom = chrom,
        ref_genome = ref_genome,
        title = title,
        highlight = loc,
        suffix = symbol,
        height = cfg$height,
        width = cfg$width
      )
    })
  }
})
