library(here)
library(tidyverse)
library(scRNAseq)
library(ensembldb)
library(scuttle)
library(cowplot)
R_SRC <- Sys.getenv("R_SRC")
source(paste0(R_SRC, "/sc_rnaseq.R"))

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
r <- qc_mads(sce, qc_spec)

sce$discard <- r$df$discard


## * Plots
# <2025-01-21 Tue> Now make the plots...
x_axis <- "cell type"
facet <- "Sex"
y_axes <- c("sum", "detected", "subsets_mito_percent")
titles <- c("Total count", "Detected features", "Percent mitochondrial reads")
library(scater)
library(ggpubr)
plots <- list()
for (i in seq_along(y_axes)) {
  y <- y_axes[i]
  title <- titles[i]
  args <- list(x = x_axis, y = y, color_by = "discard", show_violin = TRUE)
  if (!is.null(facet)) args$other_fields <- facet
  p <- do.call(\(...) plotColData(sce, ...), args) + ggtitle(title) +
    guides(color = guide_legend(title = "Discard"))
  if (!is.null(facet)) p <- p + facet_wrap(as.formula(paste("~", facet)))
  if (i != 1 && !is.null(facet)) {
    p <- p + theme(strip.text = element_blank(), strip.background = element_blank())
  }
  if (i != length(y_axes)) {
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    )
  }
  plots[[i]] <- p
}
arranged <- ggarrange(plotlist = plots, common.legend = TRUE, ncol = 1)


# Make sure that we aren't removing valid cells with both high mitochondrial percentages
# and high total counts (e.g. metabolically active cells)
mito_sum_plot <- plotColData(sce, x = "sum", "subsets_mito_percent", colour_by = "discard")

batch_check <- isOutlier(sce$sum, type = "higher", batch = sce$disease)

tb <- colData(sce_d) |> as_tibble()

tb |> ggplot(aes(
  x = factor(scDblFinder.cluster), y = cell.type,
  alpha = stat_count()
)) +
  geom_tile()
