suppressMessages({
  library(here)
  library(tidyverse)
  library(glue)
  library(ggsurvfit)
  library(AnnotationHub)
  library(ggplot2)
  library(Biobase)
  library(GEOquery)
  library(CalibrationCurves)
  library(SurvMetrics)
  library(glmnet)
  library(survival)
  library(limma)
  setAnnotationHubOption("CACHE", here(".cache", "AnnotationHub"))
  source(here("src", "R", "utils.R"))
  if (path.expand("~") == "/home/shannc") {
    reticulate::use_condaenv("stem-base")
  } else {
    reticulate::use_condaenv("nf")
  }
  set.seed(3110)
  library(reticulate)
})


GEO_CACHE <- here(".cache", "GEO")
DATE <- format(Sys.time(), "%Y-%m-%d")
CUR_DIR <- here("analyses", "brca")
ENV <- yaml::read_yaml(here(CUR_DIR, "env.yaml"))
OUTDIR <- here("analyses", "output", "brca", "re_endopredict")
dir.create(OUTDIR)
CALIB_OUTDIR <- here(OUTDIR, "calibration")
dir.create(CALIB_OUTDIR)
SC_OUTDIR <- here(OUTDIR, "survival_curves")
dir.create(SC_OUTDIR)


MPATH <- ENV$metadata_path
TIME_HORIZON <- 50

g.probe2gene <- local({
  tb <- read_tsv(here(ENV$probe_mapping)) |>
    dplyr::select(`HGNC symbol`, `AFFY HG U133A 2 probe`)
  setNames(tb$`HGNC symbol`, tb$`AFFY HG U133A 2 probe`)
})
g.gene2probe <- local({
  tb <- read_tsv(here(ENV$probe_mapping)) |>
    group_by(`HGNC symbol`) |>
    summarise(probe = list(`AFFY HG U133A 2 probe`))
  # one-to-many relationship
  setNames(tb$probe, tb$`HGNC symbol`)
})

original_set <- read_tsv(here(CUR_DIR, "endopredict_candidates.tsv")) |>
  rename_with(\(x) str_replace_all(x, " ", "_"))

# Authors already re-normalized the data from the different cohorts in the same way as theirs.
# The units are normalized such that output signal intensities of each array have a mean target intensity of 500
# shared variable is time to distant recurrence

# Unit of dmfs_time is months
shared_cols <- c("geo_accession", "subcohort", "dmfs_time", "dmfs_status")

mfiles <- list(
  GSE26971 = "GSE26971.tsv", # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26971
  GSE12093 = "GSE12093.tsv", # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12093
  `GSE6532-GPL96` = "GSE6532_LUMINAL_demo.txt", # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6532
  `GSE4922-GPL96` = "GSE4922-GPL96.tsv" # External validation set
)

expr_dir <- here(ENV$data_path, "GSE26971")
efiles <- list(
  "GSE26971_GSE12093_dataset.txt" = expr_dir,
  "GSE26971_GSE6532_dataset.txt" = expr_dir,
  "GSE26971.tsv" = expr_dir,
  "GSE4922.tsv" = here(ENV$data_path, "GSE4922-GPL96")
)


## * Gather metadata
meta <- lapply(names(mfiles), \(x) {
  file <- mfiles[[x]]
  tb <- read_tsv(here(MPATH, file))
  spec <- ENV$datasets[[x]]
  to_remap <- unlist(spec$meta_remap)
  if (!is.null(spec$id_col)) {
    tb <- dplyr::rename(tb, geo_accession = spec$id_col)
  }
  if (!is.null(to_remap)) tb <- dplyr::rename(tb, to_remap)
  if (x == "GSE6532-GPL96") {
    tb <- mutate(tb, dmfs_time = floor(dmfs_time / (365.24 / 12)))
  } else if (x == "GSE4922-GPL96") {
    tb <- mutate(tb, dmfs_time = floor(dmfs_time * 12))
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
    mutate(tb, dmfs_time = dmfs_time) |>
    dplyr::select(all_of(shared_cols))
}) |>
  bind_rows()

## * Gather expression

expr <- lmap(efiles, \(x) read_tsv(here(x[[1]], names(x[1])))) |>
  purrr::reduce(\(x, y) left_join(x, y, by = join_by("ID_REF")))
meta <- dplyr::filter(meta, geo_accession %in% colnames(expr)) |>
  column_to_rownames(var = "geo_accession")
expr <- dplyr::select(expr, all_of(c("ID_REF", rownames(meta)))) |>
  column_to_rownames(var = "ID_REF") |>
  as.matrix()
expr <- expr[, sort(colnames(expr))]
meta <- meta[sort(rownames(meta)), ]

# Represent genes with max of the probe value in group
expr <- local({
  rownames_to_column(as.data.frame(expr), var = "id") |>
    mutate(symbol = g.probe2gene[id]) |>
    dplyr::filter(!is.na(symbol)) |>
    dplyr::group_by(symbol) |>
    dplyr::summarise(across(where(is.numeric), max)) |>
    dplyr::ungroup() |>
    column_to_rownames(var = "symbol") |>
    as.matrix()
})

# ------------- Routines below only depend on expr, meta and the probe mappings

g.eset <- ExpressionSet(
  assayData = expr,
  phenoData = AnnotatedDataFrame(as.data.frame(meta))
)

g.eset <- g.eset[, !is.na(g.eset[["dmfs_status"]])]
## train_idx <- createDataPartition(g.eset[["dmfs_status"]], p = 0.8)[[1]]

g.train_eset <- g.eset[, pData(g.eset)$subcohort != "GSE4922-GPL96"]
g.test_eset <- g.eset[, pData(g.eset)$subcohort == "GSE4922-GPL96"]

## * Analyses

result <- list(
  de = list(),
  lr = list(),
  penalized_cox = list(),
  endopredict_candidates = list(
    selection = original_set$Gene_name,
    probes = original_set$Probe_set
  ),
  endopredict_rep = list()
)
other_lists <- yaml::read_yaml(here(CUR_DIR, "gene_lists.yaml"))
for (l in names(other_lists)) {
  symbols <- other_lists[[l]]
  probes <- discard(unlist(g.gene2probe[symbols]), is.na)
  result[[l]] <- list(selection = symbols, probes = probes)
}

## ** Helper functions

survival_curves <- function(model, data, cumulative = FALSE) {
  sfit <- survfit(model, data)
  key <- if (cumulative) "cumhaz" else "surv"
  sfit[[key]] |>
    as.data.frame() |>
    `rownames<-`(sfit$time)
}

bind_for_surv <- function(surv, eset) {
  if ("ExpressionSet" %in% class(eset)) {
    df <- as.data.frame(t(exprs(eset)))
  } else {
    df <- eset
  }
  cbind(as.data.frame(as.matrix(surv)), df)
}

safe_ibs <- function(model, data, set) {
  tb <- tibble(metric = "integrated_brier_score", Estimate = NA, set = set)
  try(
    tb <- tibble(
      metric = "integrated_brier_score",
      Estimate = IBS(model, data),
      set = set
    )
  )
  tb
}

##' Helper function to get scoring metrics from fitted coxph model `model`
coxph2tb <- function(model) {
  sm <- summary(model)
  tibble(
    coef = model$coefficients,
    wald_test = sm$waldtest["test"],
    wald_test_p = sm$waldtest["pvalue"],
    lr_test = sm$logtest["test"],
    lr_test_p = sm$logtest["pvalue"],
    logrank_test = sm$sctest["test"],
    logrank_test_p = sm$sctest["pvalue"]
  )
}


viz_routine <- function(embeddings, x, y, prefix) {
  with_meta <- cbind(embeddings, pData(g.eset))
  viz <- ggplot(
    with_meta,
    aes(
      x = !!as.symbol(x),
      y = !!as.symbol(y),
      color = dmfs_time,
      shape = factor(dmfs_status)
    )
  ) +
    geom_point()
  subcohort <- ggplot(
    with_meta,
    aes(x = !!as.symbol(x), y = !!as.symbol(y), color = subcohort)
  ) +
    geom_point() +
    theme(axis.title.y = element_blank())
  combined <- cowplot::plot_grid(viz, subcohort)
  ggsave(
    here(OUTDIR, glue("{prefix}_combined.png")),
    combined,
    width = 15,
    height = 8
  )
}

random_feature_performance <- function(
  x,
  y,
  feature_set,
  set_size,
  method,
  n = 5,
  trControl = NULL,
  metric = "Accuracy"
) {
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

cv_result2tb <- function(cv_result) {
  tibble(
    lambda = cv_result$lambda,
    mean_cv_error = cv_result$cvm,
    stderr = cv_result$cvsd,
  )
}

surv2structured_array <- function(surv, obj_name) {
  py_run_string("import numpy as np")
  as_mat <- as.data.frame(as.matrix(surv)) |>
    mutate(status = as.logical(status)) |>
    relocate(status, .before = time)
  py$tmp <- as_mat
  py_run_string(
    glue(
      "{obj_name} = np.array([(i['status'], i['time']) for _, i in tmp.iterrows()],
 dtype = [('s', 'bool'), ('t', 'float32')])"
    )
  )
}


## ** Visualization

## *** PCA
pca_obj <- read_existing(
  here(OUTDIR, "prcomp.rds"),
  \(f) {
    obj <- prcomp(t(exprs(g.eset)))
    saveRDS(obj, f)
    obj
  },
  readRDS
)

viz_routine(as.data.frame(pca_obj$x), "PC1", "PC2", "pca")


## *** UMAP

umap_obj <- read_existing(
  here(OUTDIR, "umap.rds"),
  \(f) {
    module <- reticulate::import("umap")
    umap <- module$UMAP()
    vals <- umap$fit_transform(t(exprs(g.eset))) |>
      `colnames<-`(c("U1", "U2")) |>
      as_tibble()
    saveRDS(vals, f)
    vals
  },
  readRDS
)

viz_routine(umap_obj, "U1", "U2", "umap")

## ** On binary outcome of no metastasis in the observed period

# NOTE: This and logistic regression shouldn't do so well

## *** DE

mm <- cbind(
  as.integer(g.train_eset[["dmfs_status"]] == 0),
  as.integer(g.train_eset[["dmfs_status"]] == 1)
) |>
  as.data.frame()
colnames(mm) <- c("no_recurrence", "recurrence")

fit <- lmFit(g.train_eset, mm)
contrast <- makeContrasts(no_recurrence - recurrence, levels = mm)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2, robust = TRUE)
top <- topTable(fit2, number = dim(g.train_eset)[1], p.value = 0.05)

result$de$selection <- rownames(top)
result$de$common <- intersect(original_set$Gene_name, top$symbol)

## *** Logistic regression with Lasso

## ** Survival analysis

g.response <- Surv(
  time = g.train_eset[["dmfs_time"]],
  event = g.train_eset[["dmfs_status"]]
)
g.test_time <- g.test_eset[["dmfs_time"]]
g.test_status <- g.test_eset[["dmfs_status"]]
g.train_time <- g.train_eset[["dmfs_time"]]
g.train_status <- g.train_eset[["dmfs_status"]]
g.test_response <- Surv(
  time = g.test_eset[["dmfs_time"]],
  event = g.test_eset[["dmfs_status"]]
)

## *** Original EndoPredict
## Method is univariate and bivariate cox

ori_wrapper <- function(eset) {
  top2a_probes <- c("201291_s_at", "201292_at")
  features <- rownames(g.train_eset)
  lapply(features, \(feature) {
    uni_data <- bind_for_surv(
      g.response,
      g.train_eset[rownames(g.train_eset) == feature, ]
    )
    uni_model <- coxph(Surv(time, status) ~ ., data = uni_data, x = TRUE)
    uni_tb <- broom::tidy(uni_model) |>
      mutate(method = "univariate")

    bi_data <- bind_for_surv(
      g.response,
      g.train_eset[rownames(g.train_eset) %in% c(feature, top2a_probes)]
    )
    bi_model <- coxph(Surv(time, status) ~ ., data = bi_data, x = TRUE)
    bi_tb <- broom::tidy(bi_model) |>
      mutate(method = "bivariate") |>
      filter(!str_remove_all(term, "`") %in% top2a_probes)
    bind_rows(uni_tb, bi_tb)
  }) |>
    bind_rows() |>
    mutate(p.adjust = p.adjust(p.value))
}

endopredict_selection <- read_existing(
  here(OUTDIR, "endopredict_original_selection.tsv"),
  \(f) {
    tb <- ori_wrapper()
    write_tsv(tb, f)
    tb
  },
  read_tsv
)

# [2025-11-17 Mon] no
ep_significant <- dplyr::filter(endopredict_selection, p.adjust <= 0.05)

result$endopredict_rep$selection <- ep_significant$term
result$endopredict_rep$common <- intersect(
  ep_significant$term,
  original_set$Gene_name
)

## *** Multivariable Cox

glmnet_cox_wrapper <- function(name, eset, penalized = FALSE) {
  read_existing(
    here(OUTDIR, glue("{name}.rds")),
    \(f) {
      x <- t(exprs(eset))
      y <- g.response
      if (penalized) {
        mod <- cv.glmnet(x, y, family = "cox", type.measure = "C", keep = TRUE)
      } else {
        mod <- bigGlm(x, y, family = "cox")
      }
      saveRDS(mod, f)
      mod
    },
    readRDS
  )
}

penalized_cox <- glmnet_cox_wrapper("penalized_cox", g.train_eset, TRUE)

coefs <- coef(penalized_cox, s = "lambda.min") |>
  as.matrix() |>
  `colnames<-`("beta") |>
  as.data.frame() |>
  rownames_to_column(var = "symbol")


nonzero <- dplyr::filter(coefs, beta > 0)
result$penalized_cox$selection <- unique(nonzero$symbol)
result$penalized_cox$probes <- nonzero$ID_REF
result$penalized_cox$common <- intersect(
  original_set$Gene_name,
  result$cox$selection
)

## *** Univariate Cox

# This is what the original endopredict authors did

## * Feature set comparison

# Compare the best C measure that other feature sets can obtain against the one
# obtained by the full penalized cox model

# NOTE: due to issues with predict.glmnet
#   (it won't work for certain models, and raises uninformative errors),
#   must resort to evaluation with scikit-survival
# TODO: SurvMetrics' calculation of IBS is super slow, can replace with something else

# Unused: if you need to go back to scikit-survival
## surv2structured_array(g.response, "y_train")
## surv2structured_array(g.test_response, "y_test")

g.times <- seq(
  min(g.train_eset[["dmfs_time"]]),
  max(g.train_eset[["dmfs_time"]]),
  by = 20
)[-1]

# %%
# Certain metrics depend on the evaluated time `timeHorizon` parameter
# - OE
# - Brier
# - Everything in `Calibration$Statistics`: ICI, E50, E90, Emax

# Helper function for getting metrics for survival analysis on train and test datasets
# Including the former gives an idea to the degree of overfitting
get_survival_metrics <- function(model, data, name) {
  cal <- valProbSurvival(
    model,
    valdata = data,
    plotCal = "ggplot",
    timeHorizon = TIME_HORIZON,
    nk = 5
  )
  with_intervals <- list(
    InTheLarge = c("Calibration", "OE"),
    TimeDependentAUC = "Uno AUC",
    Slope = c("Calibration", "calibration slope")
  )

  other_stats <- lmap(with_intervals, \(x) {
    n <- names(x)
    x <- x[[1]]
    if (length(x) == 1) {
      entry <- cal$stats[[n]]
    } else {
      entry <- cal$stats[[x[1]]][[n]]
      x <- x[2]
    }
    tibble(
      metric = x,
      Estimate = entry[x[[1]]],
      `2.5 %` = entry["2.5 %"],
      `97.5 %` = entry["97.5 %"],
    )
  }) |>
    bind_rows()

  ggsave(here(CALIB_OUTDIR, glue("{name}_curve.png")))
  bind_rows(
    other_stats,
    tibble(metric = "integrated_brier_score", Estimate = IBS(model, data)),
    rownames_to_column(as.data.frame(cal$stats$Concordance), var = "metric"),
  )
}


evaluate_cox <- function(feature_subset = NULL, name = NULL) {
  cur_eset <- g.train_eset[rownames(g.train_eset) %in% feature_subset, ]
  cur_test_eset <- g.test_eset[rownames(g.test_eset) %in% feature_subset, ]
  x_train <- bind_for_surv(g.response, cur_eset)
  x_test <- bind_for_surv(g.test_response, cur_test_eset)
  model <- coxph(Surv(time, status) ~ ., data = x_train, x = TRUE)

  ## train_metrics <- get_survival_metrics(model, x_train) |> mutate(set = "train")
  test_metrics <- get_survival_metrics(model, x_test, name) |>
    mutate(set = "test")
  metrics <- bind_rows(
    test_metrics,
    tibble(
      metric = "integrated_brier_score",
      Estimate = IBS(model, x_train),
      set = "train"
    )
  )
  tb <- mutate(metrics, name = name)
  message(glue("{name} complete"))
  return(tb)
}

## ** Overlap characteristics

overlap_df <- local({
  completed <- result[map_int(result, length) > 0]
  sets <- names(completed)
  all_symbols <- lapply(completed, \(x) x$selection) |> unlist()
  to_df <- lapply(
    all_symbols,
    \(sym) {
      sapply(completed, \(s) {
        sym %in% s$selection
      })
    }
  ) |>
    `names<-`(all_symbols)
  as.data.frame(to_df) |> mutate(across(everything(), as.integer))
})
set_dist <- vegan::vegdist(overlap_df, method = "jaccard")

## ** Run

# Include random gene sets as a baseline
with_randoms <- result
n_random_genes <- 80
n_random_iter <- 5
for (i in seq(n_random_iter)) {
  randoms <- sample(rownames(g.eset), n_random_genes)
  with_randoms[[glue("random{n_random_genes}_{i}")]] <- list(
    selection = randoms
  )
}

c_indices <- read_existing(
  here(OUTDIR, "feature_set_metrics.tsv"),
  \(f) {
    tb <- lapply(names(with_randoms), \(method) {
      if (length(with_randoms[[method]]) > 0) {
        evaluate_cox(
          feature_subset = with_randoms[[method]]$selection,
          name = method
        )
      } else {
        tibble()
      }
    }) |>
      bind_rows()
    write_tsv(tb, f)
    tb
  },
  read_tsv
)

## ** Visualize results

get_clustered_curves <- function(feature_subset, name, outdir, k = 8) {
  responses <- list(train = g.response, test = g.test_response)
  esets <- list(train = g.train_eset, test = g.test_eset)

  plots <- lapply(c("train", "test"), \(n) {
    eset <- esets[[n]]
    response <- responses[[n]]
    cur <- t(exprs(eset[rownames(eset) %in% feature_subset, ]))
    clust <- hclust(dist(cur))
    clusters <- paste0("cluster_", cutree(clust, k = k))
    obj <- bind_for_surv(response, data.frame(cluster = clusters))
    plot <- survfit2(Surv(time, status) ~ cluster, data = obj) |>
      ggsurvfit() +
      add_risktable(
        theme = list(theme(plot.title = element_text(face = "bold")))
      ) +
      ggtitle(glue("Feature set: {name}")) +
      xlab("Follow up time (months)") +
      scale_ggsurvfit() +
      ggsurvfit_build()
    if (n == "test") {
      plot <- plot + theme(axis.title.y = element_blank())
    }
    plot
  }) |>
    `names<-`(c("train", "test"))

  plot <- cowplot::plot_grid(plots$train, plots$test, ncol = 2)
  ggsave(
    filename = glue("{outdir}/{name}_curves.png"),
    plot = plot,
    height = 12,
    width = 15
  )
}

for (n in names(with_randoms)) {
  if (length(with_randoms[[n]]) > 0) {
    get_clustered_curves(with_randoms[[n]]$selection, n, outdir = SC_OUTDIR)
  }
}
