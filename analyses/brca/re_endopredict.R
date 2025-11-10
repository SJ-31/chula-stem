suppressMessages({
  library(here)
  library(tidyverse)
  library(glue)
  library(AnnotationHub)
  library(Biobase)
  library(GEOquery)
  setAnnotationHubOption("CACHE", here(".cache", "AnnotationHub"))
})
geo_cache <- here(".cache", "GEO")
date <- format(Sys.time(), "%Y-%m-%d")
cur_dir <- here("analyses", "brca")
env <- yaml::read_yaml(here(cur_dir, "env.yaml"))
outdir <- here("analyses", "output", "brca", "re_endopredict")
dir.create(outdir)
mpath <- env$metadata_path

probe2gene <- local({
  tb <- read_tsv(here(mpath, "GSE6532_LUMINAL_annot.txt")) |>
    dplyr::select(probe, HUGO.gene.symbol)
  setNames(tb$`HUGO.gene.symbol`, tb$probe)
})

original_set <- read_tsv(here(cur_dir, "endopredict_candidates.tsv")) |>
  rename_with(\(x) str_replace_all(x, " ", "_"))

# Authors already re-normalized the data from the different cohorts in the same way as theirs
# shared variable is time to distant recurrence

# Unit of dmfs_time is months
shared_cols <- c("geo_accession", "subcohort", "dmfs_time", "dmfs_status")

mfiles <- list(
  GSE26971 = "GSE26971.tsv", # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26971
  GSE12093 = "GSE12093.tsv", # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12093
  `GSE6532-GPL96` = "GSE6532_LUMINAL_demo.txt" # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6532
)

efiles <- c(
  "GSE26971_GSE12093_dataset.txt",
  "GSE26971_GSE6532_dataset.txt",
  "GSE26971.tsv"
)

expr_dir <- here(env$raw_path, "GSE26971")

## * Gather metadata
meta <- lapply(names(mfiles), \(x) {
  file <- mfiles[[x]]
  tb <- read_tsv(here(mpath, file))
  spec <- env$datasets[[x]]
  to_remap <- unlist(spec$meta_remap)
  if (!is.null(spec$id_col)) {
    tb <- dplyr::rename(tb, geo_accession = spec$id_col)
  }
  if (!is.null(to_remap)) tb <- dplyr::rename(tb, to_remap)
  if (x == "GSE6532-GPL96") {
    tb <- mutate(tb, dmfs_time = dmfs_time / (365.24 / 12))
  }
  to_replace <- spec$meta_replace
  for (col in names(to_replace)) {
    if (is.character(tb[[col]])) {
      vec <- str_to_lower(tb[[col]])
    } else {
      vec <- tb[[col]]
    }
    tb[[col]] <- unlist(to_replace[[col]])[vec] |>
      unlist(use.names = FALSE)
    if (!is.null(to_replace[[".default"]])) {
      tb[[col]] <- replace_na(tb[[col]], to_replace[[".default"]])
    }
  }
  if (!"subcohort" %in% names(to_remap)) tb <- mutate(tb, subcohort = x)
  tb |>
    mutate(tb, dmfs_time = floor(dmfs_time)) |>
    dplyr::select(all_of(shared_cols))
}) |>
  bind_rows()

## * Gather expression

expr <- lapply(efiles, \(x) read_tsv(here(expr_dir, x))) |>
  reduce(\(x, y) left_join(x, y, by = join_by("ID_REF")))
meta <- dplyr::filter(meta, geo_accession %in% colnames(expr)) |>
  column_to_rownames(var = "geo_accession")
expr <- dplyr::select(expr, all_of(c("ID_REF", rownames(meta)))) |>
  column_to_rownames(var = "ID_REF") |>
  as.matrix()
expr <- expr[, sort(colnames(expr))]
meta <- meta[sort(rownames(meta)), ]
eset <- ExpressionSet(
  assayData = expr,
  phenoData = AnnotatedDataFrame(as.data.frame(meta))
)

eset <- eset[, !is.na(eset[["dmfs_status"]])]


## * Analyses

random_feature_performance <- function(
    x,
    y,
    feature_set,
    set_size,
    method,
    n = 5,
    trControl = NULL,
    metric = "Accuracy") {
  tr <- if (is.null(trControl)) {
    trainControl(method = "cv", number = 5)
  } else {
    trControl
  }
  models <- lapply(seq(n), \(x) {
    random_set <- sample(feature_set, set_size)
    cur <- x[, colnames(x) %in% random_set]
    train(x = x, y = y, trControl = tr)
  })
  list(
    models = models,
    metric = unlist(lapply(models, \(x) x$results[[metric]][1]))
  )
}

results <- list(de = list(), lr = list(), cox = list())

## ** On binary outcome of no metastasis in the observed period

## *** DE

library(limma)
mm <- cbind(
  as.integer(eset[["dmfs_status"]] == 0),
  as.integer(eset[["dmfs_status"]] == 1)
) |>
  as.data.frame()
colnames(mm) <- c("no_recurrence", "recurrence")

fit <- lmFit(eset, mm)
contrast <- makeContrasts(no_recurrence - recurrence, levels = mm)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2, robust = TRUE)
top <- topTable(fit2, number = dim(eset)[1], p.value = 0.05)
top$symbol <- probe2gene[rownames(top)]

results$de$selection <- unique(top$symbol) |> discard(is.na)
results$de$common <- intersect(original_set$Gene_name, top$symbol)

## **** Check selection

de_filtered <- eset[rownames(expr) %in% rownames(top), ]
library(caret)

xgb_mod <- train(
  x = t(exprs(de_filtered)),
  y = factor(eset[["dmfs_status"]]),
  method = "xgbTree",
  trControl = trainControl(method = "cv", number = 5)
)

de_random_accs <- random_feature_performance(
  x = t(exprs(eset)),
  y = factor(cur[["dmfs_status"]]),
  feature_set = rownames(eset),
  method = "xgbTree",
  set_size = nrow(top),
)
# [2025-11-10 Mon] Random sets have essentially indistinguishable performance from the DE sets

## *** Logistic regression with Lasso

## ** Survival analysis

# [2025-11-10 Mon] TODO: read up on survival analysis and how to cross-validate
# you wanna see whether random gene selections will do just as well as the traditional selection process
