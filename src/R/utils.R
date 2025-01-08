read_with_filename <- function(x, col = "filename") {
  read_tsv(x) |> mutate(!!col := basename(x))
}

#' Flatten the character vector `char_vec` by `separator`, then
#'  get unique values
flatten_by <- function(char_vec, separator = ";", collapse = TRUE, unique = FALSE) {
  helper <- function(str) {
    if (is.na(str)) {
      return(NA)
    }
    lapply(str, \(x) {
      str_split_1(x, pattern = separator)
    }) |> unlist()
  }

  applied <- lapply(char_vec, helper) |> unlist()
  if (unique) {
    applied <- unique(applied)
  }
  if (collapse) {
    paste0(applied, collapse = separator)
  } else {
    list(applied)
  }
}

tb2map <- function(tb, keys, values, list = TRUE) {
  if (list) {
    as.list(tb[[values]]) |> `names<-`(tb[[keys]])
  } else {
    tb[[values]] |> `names<-`(tb[[keys]])
  }
}

into_char_list <- function(col, separator = ";") {
  lapply(col, \(x) {
    if (is.na(x)) {
      NA
    } else {
      str_split_1(x, ";")
    }
  })
}

modes <- function(x) {
  x <- discard(x, is.na)
  if (length(x) == 0) {
    return(NA)
  }
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}
