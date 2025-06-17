library(tidyverse)
library(here)
library(glue)
source(here("src", "R", "utils.R"))
source(here("analyses", "pdac", "main.R"))

mafs <- list.files(
  here(outdir, "old_maf"),
  pattern = ".*maf",
  full.names = TRUE
)


subjects <- basename(mafs) |>
  str_remove("PDAC_") |>
  str_remove("_C_filtered_variants_funco_selected.maf")

old_calls <- lapply(
  mafs,
  \(x) {
    tb <- read_tsv(x, comment = "#")
    mutate(
      tb,
      CGC_Chr = as.character(CGC_Chr),
      `HGNC_OMIM_ID(supplied_by_OMIM)` = as.character(
        `HGNC_OMIM_ID(supplied_by_OMIM)`
      ),
      gnomAD_genome_OriginalContig = as.character(gnomAD_genome_OriginalContig),
    ) |>
      mutate(across(contains("exome_AF"), as.double))
  }
) |>
  `names<-`(subjects) |>
  list_rbind(names_to = "subject")

curated_variants <- list(
  KRAS = c(
    "p.G12D",
    "p.G12C",
    "p.G12V",
    "p.G12R",
    "25209283"
  ),
  TP53 = c(
    # DBD mutations
    "p.R248Q",
    "p.R175H",
    "p.G199V",
    "p.H193R",
    "p.G199V",
    "p.H193R",
    "p.C275F",
    "p.G244D",
    "p.G245D",
    #
    "7674797"
  ),
  CDKN2A = c(
    "p.H83R",
    "p.G35R",
    "21968200"
  )
  ## SMAD4 = c( # [2025-06-17 Tue] These variants arent'present
  ##   "p.Tyr260Ter",
  ##   "51051412",
  ##   "p.Tyr430Ter"
  ## ),
  ## KMT2C = c("p.Lys2797ArgfsTer26")
)

# %%

wanted <- old_calls |>
  filter(Hugo_Symbol %in% names(curated_variants)) |>
  select(
    subject,
    Hugo_Symbol,
    Start_Position,
    End_Position,
    cDNA_Change,
    Codon_Change,
    Protein_Change,
    Variant_Classification,
    t_ref_count,
    t_alt_count,
    AF,
    DP,
  )

table_outdir <- here(outdir, "vtables")

for (sym in unique(wanted$Hugo_Symbol)) {
  filtered <- wanted |> filter(Hugo_Symbol == sym)
  if (nrow(filtered) > 0) {
    distinct(
      filtered,
      subject,
      Hugo_Symbol,
      Codon_Change,
      .keep_all = TRUE
    ) |>
      mutate(
        id = case_when(
          is.na(Protein_Change) ~ cDNA_Change,
          .default = Protein_Change
        ),
        value = paste0(t_ref_count, ",", t_alt_count, " (", round(AF, 2), ")")
      ) |>
      pivot_wider(names_from = id, id_cols = subject, values_from = value) |>
      write_tsv(here(table_outdir, glue("{sym}_old.tsv")))
  }
}
