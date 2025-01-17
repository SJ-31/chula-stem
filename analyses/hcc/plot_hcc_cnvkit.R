here::i_am("analyses/hcc/plot_hcc_cnvkit.R")
source(here::here("analyses", "hcc", "main.R"))

path <- c(
  "analyses", "data_all", "output", "HCC", "Exome", "P17",
  "variant_calling", "5-P17-Cnvkit"
)
cnvkit_dir <- do.call(here, as.list(path))

files <- list.files(cnvkit_dir, full.names = TRUE, pattern = ".*purity.*")

target_transcript <- "ENST00000397752"

read_cns <- function(x) {
  purity <- str_extract(basename(x), "(0\\.[0-9]+)_.*", 1)
  read_tsv(x) |> mutate(purity = purity)
}

all <- lapply(files, read_cns) |> bind_rows()

purity_relation <- ggplot(all, aes(x = purity, y = cn, fill = chromosome)) +
  geom_boxplot()

met_plot <- all |>
  filter(grepl(target_transcript, gene)) |>
  ggplot(aes(x = purity, y = cn)) +
  geom_bar(stat = "identity") +
  ylab("Copy number") +
  xlab("Purity") +
  labs(title = "Relationship between copy number and purity, chr7:100,951,771-142,671,152")

save_fn(met_plot, "met_purity_plot.png")
