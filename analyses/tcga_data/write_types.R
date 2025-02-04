library(here)
source(here("analyses", "tcga_data", "main.R"))

dirs <- c("TCGA_HCC", "TCGA_CHOL", "TCGA_COAD-READ")
types <- c("LIHC", "CHOL", "COADREAD")
for (i in seq_along(types)) {
  outfile <- here(outdir, glue("{str_to_lower(dirs[i])}.rds"))
  if (!file.exists(outfile)) {
    cur <- get_tcga(here(public_data, dirs[i]), types[i])
    saveRDS(cur, outfile)
  }
}
