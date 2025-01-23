## ** Cross-referenced
intogen <- read_tsv(here("analyses/data/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"))
census <- read_csv(here("analyses/data/census.csv")) |>
  rename_with(\(x) str_replace_all(x, " ", "_") |> str_remove_all("[()]")) |>
  mutate(Role_in_Cancer = map_chr(Role_in_Cancer, \(x) {
    if (is.na(x)) {
      return(x)
    }
    str_split_1(x, ",") |>
      trimws() |>
      sort() |>
      paste0(collapse = ", ")
  }))

## *** Census

with_census <- inner_join(merged, census, by = join_by(
  x$SYMBOL == y$Gene_Symbol,
  x$census_mutations == y$Mutation_Types
)) |>
  dplyr::filter((Tier == 1) &
    (Hallmark == "Yes") &
    (!is.na(Tumour_TypesSomatic)))


mutation_types <- utils$flatten_by(census$Mutation_Types, ",", collapse = FALSE) |>
  unlist() |>
  map_chr(\(x) trimws(x)) |>
  unique()
# Mis = missense
# T = translocation
# D = large deletion
# F = frameshift
# N = ???
# O = other
# A = amplification
# S = splice site
# M = mesenchymal???

# <2025-01-10 Fri> Genes shown in this plot meet the following criteria
# - have the same consequences as observed in the COSMIC census
#     from Missense, Frameshift and Deletion
# - are Tier 1, so have documented cancer activity
# - Is associated with a known somatic tumor type
# - Has a known hallmark

with_census_plot <- with_census |>
  prettify() |>
  ggplot(aes(
    x = sample, y = factor(SYMBOL, levels = names(sum_vafs)),
    fill = Role_in_Cancer, alpha = VAF
  )) |>
  vaf_heatmap() + geom_text(aes(label = round(VAF, 2))) +
  scale_y_discrete(limits = rev)

ggsave(here("analyses", "output", "pdac_cosmic_census.png"), with_census_plot, dpi = 500, width = 15)


## *** Intogen

intogen_filtered <- intogen |> filter((IS_DRIVER == "TRUE") & (ROLE != "ambiguous") &
  (TOTAL_SAMPLES >= 50))
with_intogen <- merged |>
  separate_longer_delim("Feature", ";") |>
  inner_join(intogen_filtered, by = join_by(SYMBOL, x$Feature == y$TRANSCRIPT)) |>
  group_by(SYMBOL, sample) |>
  summarise(VAF = mean(VAF))
