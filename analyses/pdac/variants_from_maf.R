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
    Genome_Change,
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

convert_protein_change <- function(protein_change) {
  # One-letter to three-letter amino acid code map
  aa_map <- c(
    A = "Ala",
    R = "Arg",
    N = "Asn",
    D = "Asp",
    C = "Cys",
    Q = "Gln",
    E = "Glu",
    G = "Gly",
    H = "His",
    I = "Ile",
    L = "Leu",
    K = "Lys",
    M = "Met",
    F = "Phe",
    P = "Pro",
    S = "Ser",
    T = "Thr",
    W = "Trp",
    Y = "Tyr",
    V = "Val",
    "*" = "Ter",
    X = "Xaa"
  )
  if (is.na(protein_change)) {
    return("NA")
  }

  # Extract parts using regular expression
  match <- regexec("^p\\.([A-Z*X])([0-9]+)([A-Z*X])$", protein_change)
  parts <- regmatches(protein_change, match)[[1]]

  if (length(parts) != 4) {
    print(protein_change)
    stop("Input must be in the form 'p.XnnnY', e.g. 'p.R71L'")
  }

  ref <- aa_map[parts[2]]
  pos <- parts[3]
  alt <- aa_map[parts[4]]

  if (is.na(ref) || is.na(alt)) {
    stop("Unknown amino acid code in input.")
  }

  paste0("p.", ref, pos, alt)
}
# %%

old <- wanted |>
  mutate(
    HGVSp = map_chr(Protein_Change, convert_protein_change),
    variant = paste0(Genome_Change, " (", HGVSp, ")"),
    symbol = Hugo_Symbol
  ) |>
  distinct(subject, variant, .keep_all = TRUE)
# Oddly there are duplicate calls for each subject...

new <- read_csv(here(outdir, "select_genes_VEP.csv")) |>
  filter(subject %in% wanted$subject) |>
  select(-INFO_TOOL_SOURCE) |>
  distinct(subject, HGVSg, .keep_all = TRUE) |> # Discard duplicates from multiple callers
  mutate(
    Genome_Change = paste0("g.", str_replace(HGVSg, ":g.", ":")),
    HGVSp = map_chr(HGVSp, \(x) {
      if (is.na(x)) {
        "NA"
      } else {
        str_remove(x, ".*:")
      }
    }),
    variant = paste0(Genome_Change, " (", HGVSp, ")"),
    found_in = case_when(variant %in% old$variant ~ "both", .default = "new"),
    symbol = SYMBOL,
  )
old$found_in <- case_when(!old$variant %in% new$variant ~ "prev")

compare_dir <- here(outdir, "compare_prev_new")

together <- full_join(
  select(old, subject, variant, found_in, symbol),
  select(new, subject, variant, found_in, symbol),
  by = join_by(subject, variant)
) |>
  mutate(
    found_in = coalesce(found_in.x, found_in.y),
    symbol = coalesce(symbol.x, symbol.y)
  )

for (sym in unique(new$symbol)) {
  filtered <- together |> filter(symbol == sym)
  filtered |>
    pivot_wider(
      names_from = subject,
      id_cols = variant,
      values_from = found_in,
      values_fill = "-"
    ) |>
    write_tsv(here(compare_dir, glue("{sym}_compare.tsv")))
  if (!all(filtered$found_in == "new")) {
    plot <- filtered |>
      ggplot(aes(x = subject, y = variant, fill = found_in)) +
      geom_tile()
    ggsave(filename = here(compare_dir, glue("{sym}_compare.png")), plot = plot)
  }
}

# %%
