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

O <- list(main = here("analyses", "output", "brca", "re_endopredict"))
O$calib <- here(O$main, "calibration")
O$sc <- here(O$main, "survival_curves")
O$assumptions <- here(O$main, "cph_assumptions")
for (outdir in O) {
  dir.create(outdir)
}


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
    tb <- mutate(tb, dmfs_time = floor(dmfs_time * 12)) |>
      dplyr::filter(er_status == "ER+")
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

g.times <- seq(
  min(g.train_eset[["dmfs_time"]]),
  max(g.train_eset[["dmfs_time"]]),
  by = 20
)[-1]


## * Analyses

result <- list(
  de = list(),
  lr = list(),
  penalized_cox = list(),
  endopredict_candidates = list(
    selection = original_set$Gene_name,
    probes = original_set$Probe_set
  ),
  endopredict_rep = list(),
  ep_rep_committee = list()
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

filter_features <- function(eset, features) {
  eset[rownames(eset) %in% features, ]
}

bind_for_surv <- function(surv, eset) {
  if ("ExpressionSet" %in% class(eset)) {
    df <- as.data.frame(t(exprs(eset)))
  } else {
    df <- eset
  }
  cbind(as.data.frame(as.matrix(surv)), df)
}

fit_coxph <- function(eset, y, subset = NULL) {
  if (!is.null(subset)) {
    eset <- eset[rownames(eset) %in% subset, ]
  }
  train <- bind_for_surv(y, eset)
  coxph(Surv(time, status) ~ ., data = train, x = TRUE)
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
    here(O$main, glue("{prefix}_combined.png")),
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

surv_summary_tibble <- function(fit, times = NULL) {
  # Compute the survival summary, optionally at specified times
  s <- summary(fit, times = times)
  n <- length(s$time)
  strata <- if (!is.null(s$strata)) {
    rep(names(s$strata), times = s$strata)
  } else {
    NULL
  }
  tb <- tibble(
    time = s$time,
    n_risk = s$n.risk,
    n_event = s$n.event,
    n_censor = s$n.censor,
    surv = s$surv,
    std_err = s$std.err,
    conf_low = s$lower,
    conf_high = s$upper
  )
  if (!is.null(strata)) {
    tb$strata <- strata
  }
  tb
}


## *** PCA
pca_obj <- read_existing(
  here(O$main, "prcomp.rds"),
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
  here(O$main, "umap.rds"),
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
g.test_response <- Surv(
  time = g.test_eset[["dmfs_time"]],
  event = g.test_eset[["dmfs_status"]]
)

train_surv <- survfit(
  Surv(time, status) ~ 1,
  data = bind_for_surv(g.response, g.train_eset)
) |>
  surv_summary_tibble(c(0, g.times))
write_tsv(train_surv, here(O$main, "train_survival.tsv"))

test_surv <- survfit(
  Surv(time, status) ~ 1,
  data = bind_for_surv(g.test_response, g.test_eset)
) |>
  surv_summary_tibble(c(0, g.times))
write_tsv(test_surv, here(O$main, "test_survival.tsv"))

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
  here(O$main, "endopredict_original_selection.tsv"),
  \(f) {
    tb <- ori_wrapper()
    write_tsv(tb, f)
    tb
  },
  read_tsv
)

ep_significant <- dplyr::filter(endopredict_selection, p.adjust <= 0.05)

result$endopredict_rep$selection <- unique(ep_significant$term)
result$endopredict_rep$common <- intersect(
  ep_significant$term,
  original_set$Gene_name
)

testph <- fit_coxph(g.train_eset, g.response, result$endopredict_rep$selection)

## *** Replicate of final algorithm

ori_final_algo <- function(
  n_members = 4,
  n_var_per_member = 2,
  wald_threshold = 0.05,
  prog_threshold = 0.05,
  tolerance = 100
) {
  library(aod)
  cohorts <- unique(pData(g.train_eset)$subcohort)
  filtered <- filter_features(g.train_eset, result$endopredict_rep$selection)
  all_features <- result$endopredict_rep$selection
  committee <- list()

  calculate_regression <- function(cur_features) {
    from_cohorts <- lapply(cohorts, \(cohort) {
      mask <- pData(filtered)$subcohort == cohort
      cur_train <- filtered[rownames(filtered) %in% cur_features, mask]
      cur_resp <- g.response[mask, ]
      cur_train <- bind_for_surv(cur_resp, cur_train)
      model <- coxph(Surv(time, status) ~ ., data = cur_train, x = TRUE)
      cf <- as.data.frame(summary(model)$coefficients)
      vcov <- vcov(model)
      wt <- map_dbl(seq_along(nrow(cf)), \(i) {
        wald.test(Sigma = vcov, b = cf$coef, Terms = i)$result$chi2[
          "P"
        ]
      })
      tibble(
        feature = rownames(cf),
        coef = cf$coef,
        var = cf$`se(coef)`,
        cohort = cohort,
        wt = p.adjust(wt)
      )
    }) |>
      bind_rows()
    sig_features <- unique(filter(from_cohorts, wt <= wald_threshold)$feature)

    lapply(cur_features, \(feat) {
      if (!feat %in% sig_features) {
        tibble(feature = feat, coef = NA, var = NA, p = 1)
      } else {
        f_tb <- dplyr::filter(from_cohorts, feature == feat)
        v <- 1 / mean(f_tb$var)
        c <- sum((f_tb$coef * (1 / f_tb$var)) * v)
        p <- 2 * pnorm(-abs(c) / sqrt(v))
        tibble(feature = feat, coef = c, var = v, p = p)
      }
    }) |>
      bind_rows()
  }

  while (length(committee) != n_members) {
    member_features <- tibble()
    rounds_unchanged <- 0
    while (nrow(member_features) != n_var_per_member) {
      if (rounds_unchanged >= tolerance) {
        stop("Failed to converge")
      }
      for (features in suppressWarnings(setdiff(
        all_features,
        member_features$feature
      ))) {
        prev_size <- nrow(member_features)

        res <- calculate_regression(features)
        best_f <- arrange(res, p) |> slice_head(n = 1)
        if (best_f$p <= wald_threshold) {
          member_features <- bind_rows(member_features, best_f)
          if (nrow(member_features) == n_var_per_member) {
            break
          }
        }

        res2 <- calculate_regression(member_features$feature)
        insig <- res2 |>
          dplyr::filter(p > prog_threshold) |>
          pluck("feature")

        member_features <- member_features |> filter(!feature %in% insig)

        cur_size <- nrow(member_features)
        if (prev_size == cur_size) {
          rounds_unchanged <- rounds_unchanged + 1
        } else {
          rounds_unchanged <- 0
        }
      }
    }
    member_res <- calculate_regression(member_features$feature)
    committee[[length(committee) + 1]] <- member_res
    all_features <- setdiff(all_features, member_res$feature)
  }
  committee
}

committee <- ori_final_algo()
result$ep_rep_committee <- list()
result$ep_rep_committee$selection <- unlist(lapply(
  committee,
  \(x) x$feature
))
result$ep_rep_committee$common <- intersect(
  original_set$Gene_name,
  result$ep_rep_committee$selection
)


## *** Multivariable Cox

glmnet_cox_wrapper <- function(name, eset, penalized = FALSE) {
  read_existing(
    here(O$main, glue("{name}.rds")),
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

# Unused: if you need to go back to scikit-survival
## surv2structured_array(g.response, "y_train")
## surv2structured_array(g.test_response, "y_test")

# %%
# Certain metrics depend on the evaluated time `timeHorizon` parameter
# - OE
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

  ggsave(here(O$calib, glue("{name}_curve.png")))
  bind_rows(
    other_stats,
    rownames_to_column(as.data.frame(cal$stats$Concordance), var = "metric"),
  )
}

get_aucs <- function(model, test, name) {
  library(survAUC)
  auc_uno <- AUC.uno(
    Surv.rsp = g.response,
    Surv.rsp.new = g.test_response,
    lpnew = predict(model, test),
    times = g.times
  )
  auc_sh <- AUC.sh(
    Surv.rsp = g.response,
    lp = predict(model),
    lpnew = predict(model, test),
    times = g.times,
    type = "cumulative"
  )
  tb <- tibble(times = g.times)
  tb[[glue("{name}_uno")]] <- auc_uno$auc
  tb[[glue("{name}_sh")]] <- auc_sh$auc
  list(at_times = tb, iauc = list(uno = auc_uno$iauc, sh = auc_sh$iauc))
}

evaluate_cox <- function(feature_subset = NULL, name = NULL) {
  cur_test_eset <- g.test_eset[rownames(g.test_eset) %in% feature_subset, ]

  model <- fit_coxph(g.train_eset, g.response, feature_subset)
  check <- cox.zph(model)
  check$table |>
    as.data.frame() |>
    rownames_to_column(var = "feature") |>
    write_tsv(here(O$assumptions, glue("{name}.tsv")))

  x_test <- bind_for_surv(g.test_response, cur_test_eset)
  x_train <- bind_for_surv(g.response, model$x)

  test_metrics <- get_survival_metrics(model, x_test, name)

  aucs <- get_aucs(model, x_test, name)

  metrics <- bind_rows(
    test_metrics,
    tibble(
      metric = c("Uno iAUC", "Song&Zhou iAUC"),
      Estimate = c(aucs$iauc$uno, aucs$iauc$sh)
    )
  )
  tb <- mutate(test_metrics, name = name)
  message(glue("{name} complete"))
  return(list(misc = tb, auc = aucs$at_times))
}


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

tmp <- read_existing(
  list(
    misc = here(O$main, "test_set_metrics.tsv"),
    auc = here(O$main, "test_set_aucs.tsv")
  ),
  \(f) {
    tibbles <- lapply(names(with_randoms), \(method) {
      if (length(with_randoms[[method]]) > 0) {
        evaluate_cox(
          feature_subset = with_randoms[[method]]$selection,
          name = method
        )
      }
    }) |>
      discard(is.null)
    result <- reduce_by_name(
      tibbles,
      fn = bind_rows,
      fn_by_name = list(auc = \(x, y) inner_join(x, y, by = join_by(times)))
    )
    write_tsv(result$misc, f$misc)
    write_tsv(result$auc, f$auc)
    result
  },
  read_tsv
)
main_metrics <- tmp$misc
aucs <- tmp$auc


## ** Visualize results

# TODO: redo this by instead stratifying patients by their hazard as obtained by the model
# nah you should just learn how to properly interpret a cox model

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
      ggtitle(n) +
      xlab("Follow up time (months)") +
      scale_ggsurvfit() |>
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
    get_clustered_curves(with_randoms[[n]]$selection, n, outdir = O$sc)
  }
}

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

## ** Compare models

yaml::write_yaml(result, here(O$main, "feature_list.yaml"))

## *** visualize previous

# TODO:

## *** compareC

#' Compare the c-indices of fitted coxph models m1, m2
#' on `data`
#'
#' @param y_true Surv object containing observed events
#' @param x1 Dataframe of shape n_samples x n_predictors
compareC_wrapper <- function(m1, m2, y_true, x1, x2 = NULL) {
  library(compareC)
  x1 <- bind_for_surv(y_true, x1)
  pred1 <- survival_curves(m1, x1)
  if (is.null(x2)) {
    x2 <- x2
  } else {
    x2 <- bind_for_surv(y_true, x2)
  }
  pred2 <- survival_curves(m2, x2)
  status <- x1$status
  times <- rownames(pred1)
  lapply(seq_len(nrow(y_true)), \(i) {
    score1 <- pred1[, i]
    score2 <- pred2[, i]
    comparison <- compareC(times, status, score1, score2)
    comparison$est.c <- NULL
    as_tibble(comparison)
  }) |>
    bind_rows() |>
    mutate(padj = p.adjust(pval)) |>
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))
}

combinations <- combn(
  names(purrr::keep_at(with_randoms, \(names) {
    map_lgl(names, \(n) {
      length(with_randoms[[n]]) > 0 &&
        (!str_detect(n, "random") || n == glue("random{n_random_genes}_1"))
    })
  })),
  2
)

cc_wrapper <- function() {
  apply(combinations, 2, \(pair) {
    x <- pair[1]
    y <- pair[2]
    f1 <- result[[x]]$selection
    f2 <- result[[y]]$selection
    m1 <- fit_coxph(g.train_eset, g.response, f1)
    m2 <- fit_coxph(g.train_eset, g.response, f2)
    tb_train <- compareC_wrapper(m1, m2, g.response, m1$x, m2$x) |>
      mutate(set = "train")
    print(tb_train)
    tb_test <- compareC_wrapper(
      m1,
      m2,
      g.test_response,
      filter_features(g.test_eset, f1),
      filter_features(g.test_eset, f2)
    ) |>
      mutate(set = "test")
    bind_rows(tb_train, tb_test) |> mutate(comparison = glue("{x} vs {y}"))
  }) |>
    bind_rows()
}

print(dim(combinations))

## BUG: this fails horribly
## concordance_comparison <- read_existing(
##   here(O$main, "concordance_comparison.tsv"),
##   \(f) {
##     tb <- cc_wrapper()
##     write_tsv(f)
##     tb
##   },
##   read_tsv
## )
