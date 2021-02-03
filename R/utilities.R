#' Read given column of tsv
#'
#' @param col_idx index of column to read
#' @param n_cols number of columns in file
#' @param tsv_file file path to .tsv file
#'
#' @return contents of the specified column in vector form
read_given_column_of_tsv <- function(col_idx, n_cols, tsv_file) {
  type_pattern <- c(rep("_", col_idx - 1), "c", rep("_", n_cols - col_idx)) %>% paste0(collapse = "")
  dplyr::pull(readr::read_tsv(file = features_fp, col_names = FALSE, col_types = type_pattern))
}

