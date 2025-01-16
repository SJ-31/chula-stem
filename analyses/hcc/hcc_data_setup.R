## * Data setup
data_dir <- here("analyses", "data_all")
id_mapping <- read_tsv(here(
  data_dir, "reference", "genomes",
  "Homo_sapiens.GRCh38.113.gene_id_mapping.tsv"
))

# Sequencer: HiSeq 2000
# <2025-01-15 Wed> It's unstranded!
tcga_data <- here(data_dir, "public_data", "TCGA_HCC")
workflow_output <- here(data_dir, "output", "HCC", "RNASEQ")

## ** TCGA data

tcga_manifest <- read_tsv(here(tcga_data, "gdc_sample_sheet.2025-01-14.tsv")) |>
  rename_with(\(x) str_replace_all(x, " ", "_"))
has_pairs <- tcga_manifest |>
  group_by(Case_ID) |>
  filter(n() >= 2) |>
  mutate(tb = map2(File_ID, File_Name, \(dir, file) {
    read_tsv(here(tcga_data, dir, file), skip = 1) |>
      select(gene_name, unstranded) |>
      rename(count = unstranded) |>
      format_counts()
  })) |>
  ungroup()

tcga_meta <- with(has_pairs, tibble(
  files = paste0(tcga_data, "/", File_ID, "/", File_Name),
  cases = has_pairs$Case_ID,
  type = Sample_Type,
  strandedness = "unstranded",
  sequencer = "HiSeq2000"
))

## ** Chula data

chula_meta <- tibble(
  cases = list.files(workflow_output, pattern = "P.*"),
  files = paste0(
    workflow_output, "/", cases, "/tumor/2-",
    cases, "_tumor-STAR_Counts.tsv"
  ),
  type = "Primary Tumor", strandedness = "reverse",
  sequencer = "NovaSeq6000"
)

chula_counts <- lapply(chula_meta$files, \(x) {
  read_tsv(x, col_names = c("gene_id", "count")) |>
    inner_join(id_mapping, by = join_by(gene_id)) |>
    select(gene_name, count) |>
    format_counts()
})

## ** Join counts

metadata <- bind_rows(tcga_meta, chula_meta) |> mutate(files = basename(files))
tb_list <- c(has_pairs$tb, chula_counts)
new_cols <- c("gene_name", metadata$files)

group_spec <- c(rep("TCGA", nrow(has_pairs)), rep("chula", length(chula_counts)))
counts <- reduce(tb_list, \(x, y) full_join(x, y, by = join_by(gene_name))) |>
  `colnames<-`(new_cols) |>
  mutate(across(where(is.numeric), \(x) replace_na(x, 0))) %>%
  DGEList(counts = ., samples = metadata, group = group_spec)

write_rds(counts, counts_file)
