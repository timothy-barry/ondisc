# These are helper functions used for testing purposes only. None of these functions are exported.

#' Create a random matrix
#'
#' All arguments optional
#'
#' @param n_row number of rows
#' @param n_col number of columns
#' @param p_zero probability an entry will be zero
#' @param matrix_values set of values from which to draw the matrix entries
#'
#' @return a randomly-generated matrix of class TsparseMatrix
create_random_matrix <- function(n_row = NULL, n_col = NULL, p_zero = 0.95, matrix_values = 1:300) {
  if (is.null(n_row)) n_row <- sample(x = 1:5000, size = 1)
  if (is.null(n_col)) n_col <- sample(x = 1:5000, size = 1)
  m <- matrix(data = sample(x = matrix_values, size = n_row * n_col, replace = TRUE), nrow = n_row, ncol = n_col)
  r <- matrix(data = stats::rbinom(n =  n_row * n_col, size = 1, prob = 1 - p_zero), nrow = n_row, ncol = n_col)
  out <- m * r
  out <- as(object = out, Class = "TsparseMatrix")
}


#' Save random matrix as a 10X object
#'
#' This function stores a matrix in .mtx format, along with features and barcodes .tsv files for the cells and genes.
#'
#' @param m a sparse Matrix object
#' @param data_dir the directory in which to store the matrix
#' @param idx (optional) an index to append to the file names
#'
#' @return the file paths to the matrix.mtx, barcodes.tsv, and features.tsv files.
save_random_matrix_as_10x <- function(m, data_dir, idx = NULL) {
  if (!dir.exists(data_dir)) dir.create(path = data_dir, recursive = TRUE)
  to_save_locs <- get_simulation_data_fps(data_dir, idx)
  # save the matrix in .mtx format.
  Matrix::writeMM(obj = m, file = to_save_locs[["mtx"]])
  # create the barcode and feature files
  cell_barcodes <- paste0("cell_", 1:ncol(m))
  gene_names <- paste0("gene_", 1:nrow(m))
  gene_ids <- paste0("ENSG000", 1:nrow(m))
  # save the files
  readr::write_tsv(x = dplyr::tibble(cell_barcodes), file = to_save_locs[["barcodes"]], col_names = FALSE)
  readr::write_tsv(x = dplyr::tibble(gene_ids, gene_names, "Gene Expression"), file = to_save_locs[["features"]], col_names = FALSE)
  # Finally, save the original R Matrix object
  saveRDS(object = m, file = to_save_locs[["r_matrix"]])
  return(to_save_locs)
}


#' Get simulation data filepaths
#'
#' Get file paths to simulation objects given an index.
#'
#' @param data_dir directory in which the simulation objects are stored
#' @param idx an index
#'
#' @return a character vector containing file paths to the simulation data
get_simulation_data_fps <- function(data_dir, idx) {
  f_names <- paste0(paste0(c("matrix", "barcodes", "features", "r_matrix", "on_disc_matrix")), if (is.null(idx)) "" else paste0("_", idx), c(".mtx", ".tsv", ".tsv", ".rds", ".rds"))
  to_save_locs <- purrr::set_names(paste0(data_dir, "/", f_names),  c("mtx", "barcodes", "features", "r_matrix", "on_disc_matrix"))
  return(to_save_locs)
}


#' Load on disc matrix and R sparse matrix
#'
#' @param data_dir simulation data directory
#' @param idx index
#'
#' @return a list containing the on_disc_matrix and the original sparse R matrix.
load_on_disc_and_mat <- function(data_dir, idx) {
  fps <- get_simulation_data_fps(data_dir, idx)
  on_disc_matrix <- fps[["on_disc_matrix"]] %>% readRDS
  r_Matrix <- fps[["r_matrix"]] %>% readRDS
  return(list(on_disc_matrix = on_disc_matrix, r_Matrix = r_Matrix))
}


#' Compare R sparse Matrix object to on_disc object on extract
#'
#' Takes a sparse R matrix, on_disc_matrix, vector of column indexes, and vector of row indexes; verifies that the row, column, and row-column subsets match.
#'
#' @param Mat an R sparse matrix object
#' @param on_disc_mat an on_disc_matrix
#' @param col_idxs column indices
#' @param row_idxs row indices
#'
#' @return NULL
compare_Mat_on_disc_extract <- function(Mat, on_disc_mat, col_idxs, row_idxs) {
  # extract sub-matrix by column
  t1 <- Mat[,col_idxs,drop=FALSE]
  t2 <- on_disc_mat[[,col_idxs]]
  testthat::expect_true(all(t1 == t2))
  # extract sub-matrix by row
  t1 <- Mat[row_idxs,,drop=FALSE]
  t2 <- on_disc_mat[[row_idxs,]]
  testthat::expect_true(all(t1 == t2))
  # extract sub-matrix by both column and row
  t1 <- Mat[row_idxs,col_idxs,drop=FALSE]
  t2 <- on_disc_mat[[row_idxs,col_idxs]]
  testthat::expect_true(all(t1 == t2))
}


#' Create synthetic data
#'
#' Generate synthetic datasets (consisting of a matrix.mtx file, a barcodes.tsv file, and a features.tsv file) and store these datasets in simulated_data_dir, with indices appended to the file names.
#'
#' @param n_datasets number of datasets to generate
#' @param simulated_data_dir directory in which to store the generated datasets
#' @param n_row number of rows in datasets (by default random)
#' @param n_col number of columns in datasets (by default random)
#' @param seed (optional) seed to set
#' @return NULL
create_synthetic_data <- function(n_datasets, simulated_data_dir, n_row = NULL, n_col = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  for (i in 1:n_datasets) {
    cat(paste0("Generating dataset ", i, ".\n"))
    m <- create_random_matrix(n_row = n_row, n_col = n_col)
    if (i %% 3 == 0) {
      n_row <- nrow(m)
      zero_row_idxs <- sample(x = 1:n_row, size = ceiling(.05 * n_row), replace = FALSE)
      m[zero_row_idxs,] <- 0
    }
    if (i %% 5 == 0) {
      n_col <- ncol(m)
      zero_col_idxs <- sample(x = 1:n_col, size = ceiling(0.05 * n_col), replace = FALSE)
      m[,zero_col_idxs] <- 0
    }
    locs <- save_random_matrix_as_10x(m, simulated_data_dir, i)
  }
  invisible()
}


#' Get test parameters
#'
#' Obtain the parameters for a run of test. The "small" test is lightweight and can easily be run on any machine on which ondisc is installed. The "large" test, by contrast, is heavier duty and takes a bit more effort to set up.
#'
#' @param test_type "small" or "big," indicating which set of tests to run.
#'
#' @return a list of parameters: synthetic_data_dir, n_datasets, n_reps_per_dataset
#' @export
get_test_parameters <- function(test_type) {
  if (test_type == "small") {
    out <- list(synthetic_data_dir = system.file("extdata", package = "ondisc"), n_datasets = 1, n_reps_per_dataset = 3, seed = 4, n_row = 300, n_col = 900)
  }
  if (test_type == "big") {
    out <- list(synthetic_data_dir = "/Users/timbarry/Box/onDisc_all/onDisc_offsite/simulated_data", n_datasets = 15, n_reps_per_dataset = 10, seed = 4, n_row = NULL, n_col = NULL)
  }
  return(out)
}


#' get test type
#'
#' @return the type of the test, either "small" or "big."
get_test_type <- function() {
  out <- c("small", "big")[1]
  return(out)
}
