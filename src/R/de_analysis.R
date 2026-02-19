library(tidyverse)
library(edgeR)

#' Create a contrast matrix for one-vs-rest comparisons
#'
#' @description
#' @return A contrast matrix of where each column is a testing contrast
#' The ith column of the matrix performs the comparison where the ith level of `objs` is
#' tested against the mean of all other levels
ovr_contrasts <- function(objs, design, prefix = "", intercept = FALSE) {
  if (!is.factor(objs)) {
    objs <- as.factor(objs)
  }
  n_levels <- nlevels(objs)
  contr <- rbind(
    matrix(1 / (1 - n_levels), n_levels, n_levels), # Average of all
    matrix(0, ncol(design) - n_levels, n_levels)
  )
  diag(contr) <- 1
  if (intercept || colnames(design)[1] == "(Intercept)") {
    contr[1, ] <- 0
  }
  rownames(contr) <- colnames(design)
  colnames(contr) <- paste0(prefix, levels(objs))
  contr
}

get_glm_fit <- function(dge, formula, id_col) {
  dge <- normLibSizes(dge)
  var_ids <- dge$genes[[id_col]]
  var_cols <- colnames(dge$genes)
  mm <- model.matrix(as.formula(formula), data = dge$samples)
  dge <- estimateDisp(dge, design = mm, robust = TRUE)
  list(fit = glmQLFit(dge, mm, robust = TRUE), mm = mm)
}


#' Helper function to run DE with edgeR
#'
#' @description
#' The form of the analysis is to compare a given level
#' of factor `group` against the mean expression of all other levels in `group`
edgeR_ovr <- function(
  dge,
  group,
  id_col = "GENEID",
  fc_cutoff = 1.5,
  p_value = 0.05,
  treat = TRUE,
  intercept = FALSE,
  extra_contrasts = NULL,
  batch_factors = NULL
) {
  if (intercept) {
    f <- paste("~", group)
  } else {
    f <- paste("~0+", group)
  }
  if (!is.null(batch_factors)) {
    f <- glue("{f}+{paste(batch_factors, collapse = ' + ')}")
  }
  tmp <- get_glm_fit(dge, f, id_col = id_col)
  fit <- tmp$fit
  mm <- tmp$fit
  ccs <- ovr_contrasts(dge$samples[[group]], mm)
  if (!is.null(extra_contrasts)) {
    ccs <- cbind(
      ccs,
      makeContrasts(contrasts = extra_contrasts, levels = colnames(mm))
    )
  }
  to_remove <- colnames(dge$genes) |> keep(\(x) x != id_col)
  compute_contrasts_separately(
    fit,
    contrasts = ccs,
    treat = treat,
    fc_cutoff = fc_cutoff,
    p_value = p_value,
    cols_remove = to_remove
  )
}


#' Helper function for calculating contrasts separately
#'
#' @param fit glmQLFit object
#' @param contrasts Matrix of contrasts
compute_contrasts_separately <- function(
  fit,
  contrasts,
  fc_cutoff,
  treat,
  p_value,
  cols_remove = c()
) {
  tests <- list()
  top <- list()
  for (i in seq_len(ncol(contrasts))) {
    name <- colnames(contrasts)[i]
    cur_contrast <- contrasts[, i]
    if (treat) {
      tests[[i]] <- glmTreat(
        fit,
        contrast = cur_contrast,
        lfc = log2(fc_cutoff)
      )
    } else {
      tests[[i]] <- glmQLFTest(fit, contrast = cur_contrast)
    }
    tests[[i]]$comparison <- name
    top[[name]] <- topTags(tests[[i]], sort.by = "PValue", p.value = p_value) |>
      as.data.frame() |>
      tibble::as_tibble() |>
      dplyr::select(-any_of(cols_remove)) |>
      dplyr::mutate(contrast = name)
  }
  num_de <- do.call(
    "cbind",
    lapply(
      lapply(tests, \(x) {
        decideTests(
          x,
          p.value = p_value,
          lfc = ifelse(treat, log2(fc_cutoff), 0)
        )
      }),
      summary
    )
  )
  list(num_de = num_de, top = bind_rows(top))
}
