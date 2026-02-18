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
  extra_contrasts = NULL
) {
  dge <- normLibSizes(dge)
  var_ids <- dge$genes[[id_col]]
  var_cols <- colnames(dge$genes)
  if (intercept) {
    f <- paste("~", group)
  } else {
    f <- paste("~0+", group)
  }
  mm <- model.matrix(as.formula(f), data = dge$samples)
  dge <- estimateDisp(dge, design = mm, robust = TRUE)
  fit <- glmQLFit(dge, mm, robust = TRUE)

  ccs <- ovr_contrasts(dge$samples[[group]], mm)
  if (!is.null(extra_contrasts)) {
    ccs <- cbind(
      ccs,
      makeContrasts(contrasts = extra_contrasts, levels = colnames(mm))
    )
  }
  to_remove <- colnames(dge$genes) |> keep(\(x) x != id_col)

  # Compute each contrast separately, as shown in the edgeR user's guide
  tests <- list()
  top <- list()
  for (i in seq_len(ncol(ccs))) {
    name <- colnames(ccs)[i]
    cur_contrast <- ccs[, i]
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
      as_tibble() |>
      select(-any_of(to_remove)) |>
      mutate(contrast = name)
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
