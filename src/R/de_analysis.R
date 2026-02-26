library(tidyverse)
library(edgeR)
library(glue)
library(checkmate)

names_if_avail <- function(x) {
  assert_atomic(x)
  names <- names(x)
  if (is.null(names)) {
    x
  } else {
    ifelse(nchar(names) == 0, x, names)
  }
}

#' Create contrasts where each grouping of levels in `groupings` is compared against
#' all other levels for the specified factor
#'
#' @param groupings A list where names are the names of contrasts
#'    and each element is a named singleton list.
#'    The name of this list is a factor name, and the values are levels to compare all
#'    other levels of the factor against
#' @param design model matrix with same name encoding used by `model.matrix`
avr_contrasts <- function(
  groupings,
  design,
  suffix = " vs Rest",
  intercept = FALSE
) {
  assert_list(groupings, min.len = 1)
  assert_matrix(design, row.names = "named", col.names = "strict")
  levels <- colnames(design)

  make_one <- function(cname, grouping) {
    assert_list(grouping, len = 1, names = "named")
    n <- length(grouping)
    factor_name <- names(grouping)[1]

    if (!any(str_starts(levels, factor_name))) {
      stop(glue(
        "The factor {factor_name} couldn't be found in the columns of design"
      ))
    }
    g1 <- grouping[[1]]
    if (is.null(cname) || nchar(cname) == 0) {
      cname <- glue("{paste0(g1, collapse = '+')} vs Rest")
    }
    pasted <- paste0(factor_name, g1)
    all_levels <- levels[str_starts(levels, factor_name)]
    if (!any(pasted %in% colnames(design))) {
      stop(glue("{pasted} not found in column names of design matrix {design}"))
    }
    avg_g1 <- 1 / length(g1)
    avg_g2 <- 1 / -(length(all_levels) - length(g1))
    contrast <- matrix(0, ncol(design), 1)
    coeffs <- dplyr::case_when(
      levels %in% pasted ~ avg_g1,
      str_starts(levels, factor_name) ~ avg_g2,
      .default = 0
    )
    contrast[, 1] <- coeffs
    rownames(contrast) <- levels
    colnames(contrast) <- cname
    if (intercept || colnames(design)[1] == "(Intercept)") {
      contrast[1, ] <- 0
    }
    contrast
  }
  lapply(
    seq_along(groupings),
    \(i) make_one(names(groupings[i]), groupings[[i]])
  ) |>
    purrr::reduce(\(x, y) cbind(x, y))
}

#' Create a contrast matrix for one-vs-rest comparisons
#'
#' @param group name of the factor that we are using for the one-vs-rest comparison
#' @param design model matrix encoded following the `model.matrix` style i.e. columns are
#' factor names concatenated with level names without space
#' @return A contrast matrix of where each column is a testing contrast
#' The ith column of the matrix performs the comparison where the ith level of
#' `objs` is tested against the mean of all other levels
#' @description
#' This only works if `group` is the baseline (first) variable or if the other variables
#' are fully encoded as a dummy variable
ovr_contrasts <- function(
  group,
  design,
  prefix = "",
  suffix = " vs Rest",
  intercept = FALSE
) {
  level_mask <- str_starts(colnames(design), group)
  n_levels <- sum(level_mask)
  if (!level_mask[1]) {
    warning(glue(
      "If variable `{group}` isn't the baseline variable, ensure that it is encoded to have all of its levels"
    ))
  }
  if (n_levels == 2) {
    warning(
      "Only two levels of factor detected, performing one vs one comparison"
    )
    one_vs_one <- TRUE
    n_levels <- 1
  } else {
    one_vs_one <- FALSE
  }
  if (!one_vs_one) {
    contrasts <- matrix(0, ncol(design), n_levels)
    contrasts[level_mask, 1:n_levels] <- 1 / (1 - n_levels)
  } else {
    contrasts <- matrix(0, ncol(design), 1)
    contrasts[level_mask, 1] <- -1
  }
  diag(contrasts) <- 1
  if (intercept || colnames(design)[1] == "(Intercept)") {
    contrasts[1, ] <- 0
  }
  rownames(contrasts) <- colnames(design)
  levels <- colnames(design)[level_mask]
  if (!one_vs_one) {
    colnames(contrasts) <- paste0(prefix, str_remove(levels, group), suffix)
  } else {
    colnames(contrasts) <- glue("{levels[1]} vs {levels[2]}")
  }
  contrasts
}


get_glm_fit <- function(dge, formula, id_col) {
  dge <- normLibSizes(dge)
  mm <- model.matrix(as.formula(formula), data = dge$samples)
  dge <- estimateDisp(dge, design = mm, robust = TRUE)
  list(fit = glmQLFit(dge, mm, robust = TRUE), mm = mm)
}

#' Helper function to run DE with edgeR
#'
#' @description
#' The form of the analysis is to compare a given level
#' of factor `group` against the mean expression of all other levels in `group`
edgeR_glm_wrapper <- function(
  dge,
  group,
  id_col = "GENEID",
  fc_cutoff = 1.5,
  p_value = 0.05,
  treat = TRUE,
  intercept = FALSE,
  avr = NULL,
  ovr = TRUE,
  contrasts = NULL,
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
  mm <- tmp$mm
  message(glue(
    'Available columns in model matrix: {paste0(colnames(mm), collapse = ", ")}'
  ))
  if (!ovr && is.null(contrasts) && is.null(avr)) {
    stop(
      "Contrasts or any-vs-rest specification must be provided if not doing one-vs-rest comparisons"
    )
  }
  if (ovr) {
    ccs <- ovr_contrasts(group = group, design = mm)
  } else {
    ccs <- NULL
  }
  if (!is.null(avr)) {
    assert_list(avr)
    ccs <- cbind(ccs, avr_contrasts(groupings = avr, design = mm))
  }
  if (!is.null(contrasts)) {
    assert_character(contrasts)
    cmatrix <- makeContrasts(contrasts = contrasts, levels = colnames(mm))
    colnames(cmatrix) <- names_if_avail(contrasts)
    ccs <- cbind(ccs, cmatrix)
  }
  to_remove <- colnames(dge$genes) |> purrr::keep(\(x) x != id_col)
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
      dplyr::select(-dplyr::any_of(cols_remove)) |>
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
  list(num_de = num_de, top = dplyr::bind_rows(top), fit = fit)
}
