library(tidyverse)
library(here)
library(glue)

refdir <- here("analyses", "data_all", "reference")

cosmic_reference <- read_csv(here(
  "analyses",
  "pdac",
  "p23_cosmic_interesting.csv"
)) |>
  distinct() |>
  filter(tier != ".") |>
  mutate(HGVSp = str_extract(HGVSp, ":(p.*)", group = 1))

observed <- read_csv(here("analyses", "pdac", "phcase_23_review.csv"))


hgvsp_recode <- function(strs, titlecase = TRUE, reverse = FALSE) {
  aa_map <- c(
    ALA = "A",
    ARG = "R",
    ASN = "N",
    ASP = "D",
    CYS = "C",
    GLN = "Q",
    GLU = "E",
    GLY = "G",
    HIS = "H",
    ILE = "I",
    LEU = "L",
    LYS = "K",
    MET = "M",
    PHE = "F",
    PRO = "P",
    SER = "S",
    THR = "T",
    TRP = "W",
    TYR = "Y",
    VAL = "V"
  )
  if (titlecase) {
    names(aa_map) <- str_to_title(names(aa_map))
  }
  if (reverse) {
    aa_map <- names(aa_map) |> `names<-`(aa_map)
  }
  capture <- paste0("(", paste0(names(aa_map), collapse = "|"), ")")
  tb <- tibble(x = strs)

  # BUG: this doesn't work perfectly, but can't think of a non-vectorized
  # workaround
  tb |>
    separate_wider_regex(
      x,
      c(
        p = "p.",
        first = "[^[:space:]0-9]+",
        loc = "[0-9]+",
        sec = "[^[:space:]0-9]+"
      ),
      too_few = "debug"
    ) |>
    print() |>
    mutate(
      first = ifelse(str_detect(first, capture), aa_map[first], first),
      sec = ifelse(str_detect(sec, capture), aa_map[sec], sec)
    ) |>
    unite("combined", p, first, loc, sec, x_remainder, sep = "") |>
    pluck("combined")
}

if ("change_recoded" %notin% colnames(observed)) {
  observed$change_recoded <- hgvsp_recode(observed$change, reverse = TRUE)
  write_csv(observed, here("analyses", "pdac", "phcase_23_review.csv"))
}

joined <- observed |>
  inner_join(
    cosmic_reference,
    by = join_by(x$change_recoded == y$HGVSp, x$symbol == y$gene)
  ) |>
  distinct(id, change_recoded, .keep_all = TRUE) |>
  filter(vaf > 10)

write_csv(joined, here("analyses", "pdac", "phcase_23_confirmed_cosmic.csv"))

cosmic_classification <- read_tsv(here(
  refdir,
  "variants/Cosmic_Classification_Tsv_v104_GRCh38.tar"
))

rel_genes <- yaml::read_yaml(here(
  "analyses",
  "pdac",
  "snakemake",
  "vaf_config.yaml"
))$wanted_genes |>
  names()

joined |> filter(symbol %in% rel_genes)
