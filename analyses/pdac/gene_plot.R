library(ensembldb)
library(here)
library(reticulate)
library(paletteer)
use_condaenv("base")
source(here::here("analyses", "pdac", "main.R"))

db <- EnsDb(here("analyses", "data", "Homo_sapiens.GRCh38.113.sqlite"))
seqlevelsStyle(db) <- "UCSC"

vep_files <- list.files(
	data_path,
	pattern = "8-P[0-9_]+-VEP_small.tsv$",
	recursive = TRUE,
	full.names = TRUE
)

names <- utils$basename_no_ext(vep_files) |>
	map_chr(\(x) str_extract(x, "8-(P[0-9_]+)-VEP_small", 1))

P <- new.env()
source(here("src", "R", "plotting.R"), local = P)

P$plot_sample_variants(
	db,
	vep_files,
	"KRAS",
	here(outdir, "kras_sample_plot.png"),
	canonical_tx = NULL,
	sample_names = names,
	palette = "vapoRwave::vapoRwave"
)

P$plot_sample_variants(
	db,
	vep_files,
	"TP53",
	here(outdir, "tp53_sample_plot.png"),
	canonical_tx = NULL,
	sample_names = names,
	palette = "vapoRwave::vapoRwave"
)

P$plot_sample_variants(
	db,
	vep_files,
	"CDKN2A",
	here(outdir, "cdk2na_sample_plot.png"),
	canonical_tx = NULL,
	sample_names = names,
	palette = "vapoRwave::vapoRwave"
)
