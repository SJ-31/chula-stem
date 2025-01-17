#!/usr/bin/env Rscript

if (sys.nframe() == 0) {
  library("dupRadar")
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-b", "--input_bam"),
    type = "character",
    help = "Input bam file with marked duplicates"
  )
  parser <- add_option(parser, c("-g", "--gtf"), type = "character", help = "Reference gtf")
  parser <- add_option(parser, c("-o", "--output"),
    type = "character", help = "Output file name",
    default = "dupradar.tsv"
  )
  parser <- add_option(parser, c("-t", "--threads"), type = "character", help = "N threads to use", default = 1)
  parser <- add_option(parser, c("-s", "--strandedness"),
    type = "character", help = "Strandedness - forward|reverse|unstranded"
  )
  parser <- add_option(parser, c("-u", "--unpaired"),
    type = "logical", help = "Reads were unpaired", action = "store_true"
  )
  parser <- add_option(parser, "--box_plot",
    type = "character",
    help = "Path to box plot (not drawn if left empty)", default = NULL
  )
  parser <- add_option(parser, "--density_plot",
    type = "character",
    help = "Path to density plot (not drawn if left empty)", default = NULL
  )
  args <- parse_args(parser)
  paired <- !args$unpaired
  dm <- dupRadar::analyzeDuprates(
    args$bam, args$gtf, args$strandedness, paired, args$threads
  )
  write.table(dm, args$output, sep = "\t", quote = FALSE, row.names = FALSE)
  if (!is.null(args$box_plot)) {
    png(filename = args$box_plot, quality = 100)
    duprateExpBoxplot(dm)
    dev.off()
  }
  if (!is.null(args$density_plot)) {
    png(filename = args$density_plot, quality = 100)
    duprateExpDensPlot(dm)
    dev.off()
  }
}
