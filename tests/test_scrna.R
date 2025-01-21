library(here)
library(tidyverse)
library(scRNAseq)
library(ensembldb)
library(scuttle)

here::i_am("tests/test_scrna.R")
AnnotationHub::setAnnotationHubOption("CACHE", here(".cache", "AnnotationHub"))
ah <- AnnotationHub::AnnotationHub()
data <- scRNAseq::LawlorPancreasData()
refdir <- here("analyses", "data_all", "reference")

gff_file <- here(refdir, "genomes", "Homo_sapiens.GRCh38.113.sqlite")
db <- ensembldb::EnsDb(gff_file)

anno <- AnnotationDbi::select(db,
  keys = rownames(data), keytype = "GENEID",
  columns = c("SYMBOL", "GENEBIOTYPE", "SEQNAME")
)
rowData(data) <- left_join(data.frame(GENEID = rownames(data)), anno, by = join_by(GENEID))
is_mito <- (rowData(data)$SEQNAME == "MT") %>% replace_na(FALSE)
sce <- scater::addPerCellQC(data, subsets = list(mito = is_mito))

# <2025-01-21 Tue> Are there marker genes that when expressed are indicative of poor
# cell quality? Would be easy to implement

# Do QC
qc_spec <- list(sum = list(type = "both", nmads = 5), detected = list(type = "high", nmads = 2))
r <- qc_wrapper(sce, qc_spec)
sce$discard <- r$df$discard

# <2025-01-21 Tue> Now make the plots...
x_axis <- "cell type"
facet <- "Sex"
y_axes <- c("sum", "detected", "subsets_mito_percent")
titles <- c("Total count", "Detected features", "Percent mitochondrial reads")
library(scater)
plots <- list()
for (i in seq_along(y_axes)) {
  y <- y_axes[i]
  title <- titles[i]
  args <- list(x = x_axis, y = y, color_by = "discard")
  if (!is.null(facet)) args$other_fields <- facet
  p <- do.call(\(...) plotColData(sce, ...), args) + ggtitle(title)
  plots[[i]] <- p
}

sum_plot <- plotColData(sce, x = x_axis, "sum", colour_by = "discard")
detected_plot <- plotColData(sce, x = x_axis, "detected", colour_by = "discard")
mito_plot <- plotColData(sce, x = x_axis, "subsets_mito_percent", colour_by = "discard")


if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-m", "--n_mads"),
    default = 3, help = "Minimum number of MADs to consider a cell as an outlier"
  )
  parser <- add_option(parser, c("-d", "--direction"), type = "character", help = "Directionality of MADs to consider")
  parser <- add_option(parser, c("-i", "--input"), type = "character", help = "Input counts")
  parser <- add_option(parser, c("-o", "--output"), type = "character", help = "Output file name")
  args <- parse_args(parser)
}
