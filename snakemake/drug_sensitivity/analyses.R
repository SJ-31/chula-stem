library(tidyverse)
library(glue)

## * Utils

#' Try to evaluate `expr`, returning a list of objects from
#' the attempt
#'
try_expr <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      err <<- e
      NULL
    }),
    warning = function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warning = warn, error = err)
}

#' Group samples according to specified response quantiles
#'
add_response_group <- function(
    tb,
    response_col,
    group_spec = list(
      sensitive = "<0.25",
      resistant = ">0.75",
      intermediate = NULL
    ),
    col_added = NULL) {
  if (is.null(col_added)) {
    col_added <- glue("{response_col}_group")
  }
  response <- tb[[response_col]]
  cdf <- ecdf(response)
  make_na <- is.na(response)
  to_check <- discard(group_spec, is.null)
  fill_group <- names(keep(group_spec, is.null))[1]

  checked <- lapply(names(to_check), \(g) {
    result <- list()
    spec <- group_spec[[g]]
    val <- cdf(response)
    vec <- eval(parse(text = glue("val{group_spec[[g]]}")))
    result[[g]] <- ifelse(vec, g, NA)
    as_tibble(result)
  }) |>
    bind_cols()
  new_col <- coalesce(!!!checked)
  new_col <- ifelse(is.na(new_col) & !make_na, fill_group, new_col)

  tb[[col_added]] <- new_col
  tb
}

## * Functions

#' For each variable (columns of `tb`), perform pairwise tests to
#'  see if there is a statistically significant difference in `response`
#'  between the subgroup of `tb` with and without the feature
#'
#' @param tb tibble where samples are rows and ALL features are columns
#' @param response response vector, should be aligned to samples in `tb`
#' @param test hypothesis testing function with two arguments. First argument is
#'      the response vector for the feature `present`
#' @param correction multiple-testing correction function, takes vector of p-values as argument
binary_feature_analysis <- function(
    tb,
    response,
    test = \(x, y) wilcox.test(x, y, alternative = "greater"), # We are interested that the presence of a mutation
    correction = p.adjust) {
  if (is.character(response) && length(response) == 1) {
    tmp <- response
    response <- tb[[response]]
    tb <- select(tb, -all_of(tmp))
  }
  if (length(response) != dim(tb)[1]) {
    stop("`response` vector and sample `tb` aren't aligned!")
  }
  lapply(colnames(tb), \(feature) {
    mask <- tb[[feature]] > 0
    present <- response[mask]
    absent <- response[!mask]
    result <- list(feature = feature)
    test <- try_expr(test(present, absent))
    if (is.null(test$value)) {
      result$p_value <- NA
      result$success <- FALSE
      result$reason <- test$error$message
    } else {
      result$p_value <- test$value$p.value
      result$success <- TRUE
    }
    as_tibble(result)
  }) |>
    bind_rows() |>
    mutate(p_adjust = correction(p_value)) |>
    relocate(p_adjust, .after = p_value)
}
