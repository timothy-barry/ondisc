#' Get mtx metadata
#'
#' @param mtx_fp filepath to the mtx file
#' @param n_rows_with_comments number of rows with comments (at top of file)
#'
#' @return a list containing (i) n_genes, (ii) n_cells, (iii) the numer of data points (i.e., fraction of entries that are zero), (iv) (TRUE/FALSE) matrix is logical
get_mtx_metadata <- function(mtx_fp, n_rows_with_comments) {
  metadata <- utils::read.table(file = mtx_fp, nrows = 1, skip = n_rows_with_comments, header = FALSE, sep = " ", colClasses = c("integer", "integer", "integer"))
  n_features <- metadata %>% dplyr::pull(1)
  n_cells <- metadata %>% dplyr::pull(2)
  n_data_points <- metadata %>% dplyr::pull(3)
  if (n_data_points > .Machine$integer.max) stop("Numer of rows exceeds maximum value.")
  first_row <- utils::read.table(file = mtx_fp, nrows = 1, skip = n_rows_with_comments + 1, header = FALSE, sep = " ", colClasses = "integer")
  is_logical <- ncol(first_row) == 2
  return(list(n_features = n_features, n_cells = n_cells, n_data_points = n_data_points, is_logical = is_logical))
}


#' Get metadata for features.tsv file
#'
#' Gets metadata from a features.tsv file. As a side-effect, if MT genes are present, puts into the bag_of_variables a logical vector indicating the positions of those genes.
#'
#' @param features_fp file path to features.tsv file
#' @param bag_of_variables the bag of variables to which to add the mt_genes logical vector (if applicable)
#'
#' @return a list containing elements feature_names (logical), n_cols (integer), and wheter MT genes are present (logical)
get_features_metadata <- function(features_fp, bag_of_variables) {
  first_row <- readr::read_tsv(file = features_fp, n_max = 1, col_names = FALSE, col_types = readr::cols())
  n_cols <- ncol(first_row)
  feature_names <- ncol(first_row) >= 2
  mt_genes_present <- FALSE
  if (feature_names) {
    gene_names <- read_given_column_of_tsv(col_idx = 2, n_cols = n_cols, tsv_file = features_fp)
    mt_genes <- grepl(pattern = "^MT-", x = gene_names)
    if (any(mt_genes)) {
      mt_genes_present <- TRUE
      bag_of_variables[[arguments_enum()$mt_gene_bool]] <- mt_genes
    }
  }
  return(list(feature_names = feature_names, n_cols = n_cols, mt_genes_present = mt_genes_present))
}


#' Generate on disc_matrix_name
#'
#' Generates the name of an on_disc_matrix object given a directory. This function searches for files named on_disc_matrix_x.h5 in the specified directory. If none exists, it returns on_disc_matrix_1.h5. Else, it returns n_disc_matrix_x.h5 with a unique integer in place of x.
#'
#' @param on_disc_dir directory in which to store the on_disc_matrix.
#' @return a new name for an on_disc_matrix.
generate_on_disc_matrix_name <- function(on_disc_dir) {
  fs <- list.files(on_disc_dir)
  base_name <- "ondisc_matrix_"
  idxs <- grep(pattern = paste0(base_name, "[0-9]+.h5"), x = fs)
  if (length(idxs) == 0) {
    name <- paste0(base_name, "1.h5")
  } else {
    existing_names <- fs[idxs]
    ints_in_use <- gsub(pattern = paste0(base_name, "(\\d+).h5"), replacement = "\\1", x = existing_names) %>% as.integer()
    new_int <- max(ints_in_use) + 1
    name <- paste0(base_name, new_int, ".h5")
  }
  return(name)
}


#' Read given column of tsv
#'
#' @param col_idx index of column to read
#' @param n_cols number of columns in file
#' @param tsv_file file path to .tsv file
#'
#' @return contents of the specified column in vector form
read_given_column_of_tsv <- function(col_idx, n_cols, tsv_file) {
  type_pattern <- c(rep("_", col_idx - 1), "c", rep("_", n_cols - col_idx)) %>% paste0(collapse = "")
  dplyr::pull(readr::read_tsv(file = tsv_file, col_names = FALSE, col_types = type_pattern))
}


#' get n rows with comments mtx
#'
#' @param mtx_fp file path to mtx
#'
#' @return the number of rows in the mtx with comments
get_n_rows_with_comments_mtx <- function(mtx_fp) {
  n_rows_with_comments <- 0
  repeat {
    curr_row <- utils::read.table(mtx_fp, nrows = 1, skip = n_rows_with_comments, header = FALSE, sep = "\n") %>% dplyr::pull()
    is_comment <- substr(curr_row, start = 1, stop = 1) == "%"
    if (!is_comment) {
      break()
    } else {
      n_rows_with_comments <- n_rows_with_comments + 1
    }
  }
  return(n_rows_with_comments)
}
