read_with_filename <- function(x, col = "filename") {
  read_tsv(x) |> mutate(!!col := basename(x))
}

#' Flatten the character vector `char_vec` by `separator`, then
#'  get unique values
flatten_by <- function(char_vec, separator = ";", collapse = TRUE) {
  helper <- function(str) {
    if (is.na(str)) {
      return(NA)
    }
    lapply(str, \(x) {
      str_split_1(x, pattern = separator)
    }) |> unlist()
  }
  applied <- lapply(char_vec, helper) |>
    unlist() |>
    unique()
  if (collapse) {
    paste0(applied, collapse = separator)
  } else {
    list(applied)
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
