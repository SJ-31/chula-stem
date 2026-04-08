library(tidyverse)
library(jsonlite)
library(here)

datadir <- here("analyses", "data", "mouse_STR")
workdir <- here("snakemake", "mouse_authentication")
tb <- read_csv(here(datadir, "mouse_STR_multiplex_pcr_primers.csv"))
catalog_file <- here(workdir, "mm39_mouse_str_catalog.json")


# Repeat motifs for D4S408, D8S1106 were obtained from Table 3. of
# from https://pmc.ncbi.nlm.nih.gov/articles/PMC3221628/
repeat_tb <- read_csv(here(datadir, "mouse_STR_loci_info.csv")) |>
  add_row(
    STR_Marker = "D4S2408",
    Repeat_Motif = "[ATCT]nACCC[ATCT]n[ACCT]nATCT"
  ) |>
  add_row(
    STR_Marker = "D8S1106",
    Repeat_Motif = "[ATAG]n"
  )

bed <- tb |>
  mutate(
    chr = str_extract(
      `Chromosome Location`,
      pattern = "(NC_[0-9.]+):",
      group = 1
    ),
    bounds = str_extract(`Chromosome Location`, "[0-9]+-[0-9]+"),
  ) |>
  separate_wider_delim(bounds, "-", names = c("start", "end")) |>
  select(chr, start, end, Marker) |>
  mutate(end = end + 550) |> # To accomodate for size increase due to the repeats
  inner_join(
    select(repeat_tb, STR_Marker, Repeat_Motif),
    by = join_by(x$Marker == y$STR_Marker)
  )

to_expansion_hunter <- bed |>
  apply(1, \(row) {
    ref_region <- paste0(row["chr"], ":", row["start"], "-", row["end"])
    list(
      LocusId = unbox(row["Marker"]),
      ReferenceRegion = unbox(ref_region),
      VariantType = unbox("Repeat"),
      LocusStructure = unbox(
        str_replace_all(row["Repeat_Motif"], "\\[", "(") |>
          str_replace_all("\\]", ")") |>
          str_replace_all("n", "*")
      )
    )
  })

to_expansion_hunter |>
  write_json(catalog_file, pretty = TRUE, autounbox = TRUE)
