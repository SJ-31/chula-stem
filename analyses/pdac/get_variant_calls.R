library(tidyverse)
library(here)
library(glue)
source(here("analyses", "pdac", "main.R"))
vep_vcfs <- list.files(
  data_path,
  pattern = "7-P[0-9_]+-VEP_small.vcf.gz$",
  recursive = TRUE,
  full.names = TRUE
)
mutect_vcfs <- list.files(
  data_path,
  pattern = "5-P[0-9_]+-Mutect2_filtered.vcf.gz",
  recursive = TRUE,
  full.names = TRUE
)

bams <- list.files(
  data_path,
  pattern = "4-P[0-9_]+_tumor-recal.bam$",
  recursive = TRUE,
  full.names = TRUE
)

wanted <- c(
  KRAS = "12:25205246-25250936",
  TP53 = "17:7661779-7687546",
  MUC5B = "11:1223066-1262172",
  KMT2C = "7:152134922-152436644",
  ARID1A = "1:26693236-26782104",
  SMAD4 = "18:51028528-51085045",
  GLI3 = "7:41960949-42264100",
  CDKN2A = "9:21967752-21995301"
)
target_file <- here("analyses", "pdac", "target_gene_regions.txt")
query_str <- "[%AD\t%AF\t%GT\t%PS]\t%INFO/DP\t%INFO/SOURCE"

get_wanted_genes_mutect <- function(file) {
  args <- c("view", glue("--targets-file {target_file}"), file)
  pref <- utils$basename_no_ext(file) |> str_remove(".vcf")
  outfile <- glue("{pref}.tsv")
  cur_outdir <- here(outdir, "mutect_tsv")
  filtered <- here(cur_outdir, outfile)
  mutect_calls <- system2("bcftools", args, stdout = TRUE)
  system2(
    "bcftools",
    args = c("query", glue("-f '%CHROM\t%POS\t%REF\t%ALT\t{query_str}'")),
    input = mutect_calls,
    stdout = here(cur_outdir, "temp.tsv")
  )
  read_tsv(
    here(cur_outdir, "temp.tsv"),
    col_names = c(
      "CHROM",
      "POS",
      "REF",
      "ALT",
      "AD",
      "AF",
      "GT",
      "PS",
      "INFO_DP",
      "SOURCE"
    )
  ) |>
    select(-SOURCE) |>
    write_tsv(filtered)
}

get_wanted_genes <- function(file, vep = TRUE) {
  args <- c("view", glue("--targets-file {target_file}"), file)
  pref <- utils$basename_no_ext(file) |> str_remove(".vcf")
  outfile <- glue("{pref}.vcf")
  filtered <- here(outdir, "vcfs", outfile)
  system2("bcftools", args, stdout = filtered)

  final_file <- here(outdir, "vcfs", glue("{pref}.tsv"))
  args2 <- c(
    glue("-i {filtered}"),
    glue("-o {final_file}"),
    glue("-q '{query_str}'")
  )
  if (vep) {
    system2(here("src", "bash", "get_vep_anno.bash"), args = args2)
  }
}

get_wanted_bam <- function(file) {
  args <- c("view", "-b", "-h", glue("-L {target_file}"), file)
  pref <- utils$basename_no_ext(file)
  bam_out <- here(outdir, "bams", glue("{pref}.bam"))
  system2("samtools", args, stdout = bam_out)
}

if (path.expand("~") != "/home/shannc") {
  ## tmp <- lapply(vep_vcfs, get_wanted_genes)
  ## tmp <- lapply(mutect_vcfs, get_wanted_genes_mutect)
  ## tmp <- lapply(bams, get_wanted_bam)
}

get_mutect_inconsistent <- function(file) {
  cur_outdir <- here(outdir, "mutect_tsv")
  files <- list.files(cur_outdir)
  subjects <- str_extract(files, "5-(P[0-9_]+)-.*", group = 1)
  lapply(files, \(x) {
    read_tsv(here(cur_outdir, x)) |> mutate(AF = as.double(AF))
  }) |>
    `names<-`(subjects) |>
    list_rbind(names_to = "subject") |>
    filter(grepl("^0,", AD) & AF != 1) |>
    write_tsv(file)
}

mutect_inconsistent <- read_existing(
  here(
    outdir,
    "mutect2_inconsistent_gt_ad.tsv"
  ),
  get_mutect_inconsistent,
  read_tsv
)


all_vcfs <- list.files(
  here(outdir, "vcfs"),
  pattern = "*.tsv",
  full.names = TRUE
)
subjects <- utils$basename_no_ext(all_vcfs) |>
  str_extract("-(P.*)-VEP_small.*", group = 1)

cosmic_all <- read_tsv(here(outdir, "cosmic_metadata.tsv"))

cosmic_md <- read_tsv(here(outdir, "cosmic_metadata.tsv")) |>
  mutate(PUBMED_PMID = as.character(PUBMED_PMID)) |>
  filter(MUTATION_SOMATIC_STATUS != "Variant of unknown origin") |>
  mutate(
    study = case_when(
      is.na(COSMIC_STUDY_ID) ~ PUBMED_PMID,
      .default = COSMIC_STUDY_ID
    )
  ) |>
  group_by(HGVSG) |>
  summarise(
    n_studies = length(unique(study)),
    n_samples = length(unique(COSMIC_SAMPLE_ID))
  )

combined_anno <- lapply(all_vcfs, \(x) {
  read_tsv(x) |>
    mutate(
      Protein_position = as.character(Protein_position),
      AF = as.double(AF),
      PS = as.double(PS),
      INFO_DP = as.double(INFO_DP)
    )
}) |>
  `names<-`(subjects) |>
  list_rbind(names_to = "subject") |>
  left_join(cosmic_md, by = join_by(x$HGVSg == y$HGVSG))


combined_anno |> write_csv(here(outdir, "select_genes_VEP.csv"))

for (n in names(wanted)) {
  combined_anno |>
    filter(SYMBOL == n) |>
    group_by(subject, HGVSg) |>
    rename(n_studies_cosmic = n_studies) |>
    select(
      subject,
      HGVSg,
      HGVSp,
      n_studies_cosmic,
      Protein_position,
      CLIN_SIG,
      Consequence,
      GT,
      AD,
      AF,
      INFO_DP,
      INFO_SOURCE,
      Amino_acids,
      Existing_variation
    ) |>
    arrange(subject, HGVSg) |>
    rename(Codon_number = Protein_position) |>
    distinct() |>
    write_tsv(here(outdir, "formatted_tsv", glue("{n}.tsv")))
}

# TODO: go through the literature Aj gives to you and define the variant list here
# then recreate the plot using those variants

# TODO: suggest if you should highlight the different codons affected, since
# for some genes, each patient has a different set of protein-coding variants

combined_anno |>
  distinct() |>
  filter(!is.na(Protein_position)) |>
  ggplot(aes(x = subject, y = SYMBOL)) +
  geom_tile()


cur_symbol <- "KRAS"
current <- combined_anno %>%

  filter(SYMBOL == cur_symbol) %>%
  distinct(subject, HGVSc, .keep_all = TRUE) |>
  group_by(HGVSc) |>
  summarize(
    across(
      all_of(c(
        "CLIN_SIG",
        "n_studies",
        "HGVSg",
        "HGVSp",
        "Consequence",
        "Existing_variation"
      )),
      first
    ),
    n_subjects = n(),
  ) |>
  mutate(
    COSMIC = lapply(Existing_variation, \(x) {
      if (is.na(x)) {
        return(x)
      }
      str_split_1(x, "&") |>
        discard(\(str) str_detect(str, "^rs")) |>
        paste0(collapse = ",")
    }) |>
      unlist()
  ) |>
  arrange(desc(n_subjects))
