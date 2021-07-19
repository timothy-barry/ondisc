#' Create synthetic data
#'
#' Generate synthetic datasets (consisting of a matrix.mtx file, a barcodes.tsv file, and a features.tsv file) and store these datasets in simulated_data_dir, with indices appended to the file names.
#'
#' @param n_datasets number of datasets to generate
#' @param simulated_data_dir directory in which to store the generated datasets
#' @param n_row number of rows in datasets (by default random)
#' @param n_col number of columns in datasets (by default random)
#' @param seed (optional) seed to set
#' @param idx_start index at which to start (default 1)
#' @param logical_mat_p probability a matrix is logical
#' @return NULL
create_synthetic_data <- function(n_datasets, simulated_data_dir, n_row = NULL, n_col = NULL, seed = NULL, idx_start = 1L, logical_mat_p = 0.5) {
  if (!is.null(seed)) set.seed(seed)
  out <- vector(mode = "list", length = n_datasets)
  for (i in seq(idx_start, idx_start + n_datasets - 1L)) {
    if (n_datasets > 1) cat(paste0("Generating dataset ", i, ".\n"))
    if (i == 1) {
      m <- create_random_matrix(n_row = n_row, n_col = n_col, logical_mat = TRUE)
    } else if (i == 2) {
      m <- create_random_matrix(n_row = n_row, n_col = n_col, matrix_values = 1:10)
    } else if (stats::rbinom(1, 1, logical_mat_p)) {
      m <- create_random_matrix(n_row = n_row, n_col = n_col, logical_mat = TRUE)
    } else {
      m <- create_random_matrix(n_row = n_row, n_col = n_col, matrix_values = 1:10)
    }
    if (stats::rbinom(1, 1, 0.5)) {
      n_row <- nrow(m)
      zero_row_idxs <- sample(x = seq(1, n_row), size = ceiling(.25 * n_row), replace = FALSE)
      m[zero_row_idxs,] <- if (is.logical(m@x[1])) FALSE else 0
    }
    if (stats::rbinom(1, 1, 0.5)) {
      n_col <- ncol(m)
      zero_col_idxs <- sample(x = 1:n_col, size = ceiling(0.05 * n_col), replace = FALSE)
      m[,zero_col_idxs] <- if (is.logical(m@x[1])) FALSE else 0
    }
    # create the barcode and feature files
    cell_barcodes <- paste0("cell_", 1:ncol(m))
    gene_names <- paste0("gene_", 1:nrow(m))
    # set 1/10 of entries to MT-*
    idxs <- sort(sample(x = 1:nrow(m), size = floor(nrow(m)/10), replace = FALSE))
    gene_names[idxs] <- paste0("MT-", idxs)
    gene_ids <- paste0("ENSG000", 1:nrow(m))
    features_df <- data.frame(gene_ids, gene_names)
    # write to disk
    save_random_matrix_as_10x(m = m, data_dir = simulated_data_dir,
                              cell_barcodes = cell_barcodes, gene_names = gene_names, gene_ids = gene_ids,
                              idx = i, save_r_matrix = FALSE)
    # prepare output
    out[[i - idx_start + 1L]] <- list(r_matrix = m, cell_barcodes = cell_barcodes, features_df = features_df)
  }
  return(out)
}


#' Create a random matrix
#'
#' All arguments optional
#'
#' @param n_row number of rows
#' @param n_col number of columns
#' @param p_zero probability an entry will be zero
#' @param matrix_values set of values from which to draw the matrix entries
#' @param logical_mat should the matrix be logical (as opposed to numeric)?
#'
#' @return a randomly-generated matrix of class TsparseMatrix
#' @noRd
create_random_matrix <- function(n_row = NULL, n_col = NULL, p_zero = 0.95, matrix_values = 1:10, logical_mat = FALSE) {
  if (is.null(n_row)) n_row <- sample(x = 200:1000, size = 1)
  if (is.null(n_col)) n_col <- sample(x = 200:1000, size = 1)
  r <- matrix(data = stats::rbinom(n =  n_row * n_col, size = 1, prob = 1 - p_zero), nrow = n_row, ncol = n_col)
  if (!logical_mat) {
    m <- matrix(data = sample(x = matrix_values, size = n_row * n_col, replace = TRUE), nrow = n_row, ncol = n_col)
    out <- m * r
  } else {
    out <- r == 1
  }
  return(Matrix::Matrix(data = out, sparse = TRUE))
}


#' Load on disc matrix and R sparse matrix
#'
#' @param data_dir simulation data directory
#' @param idx index
#'
#' @return a list containing the on_disc_matrix and the original sparse R matrix.
#' @noRd
load_on_disc_and_mat <- function(data_dir, idx) {
  fps <- get_simulation_data_fps(data_dir, idx)
  on_disc_matrix <- fps[["on_disc_matrix"]] %>% readRDS
  r_Matrix <- fps[["r_matrix"]] %>% readRDS
  return(list(on_disc_matrix = on_disc_matrix, r_Matrix = r_Matrix))
}


#' Get simulation data filepaths
#'
#' Get file paths to simulation objects given an index.
#'
#' @param data_dir directory in which the simulation objects are stored
#' @param idx an index
#'
#' @return a character vector containing file paths to the simulation data
#' @noRd
get_simulation_data_fps <- function(data_dir, idx) {
  f_names <- paste0(paste0(c("matrix", "barcodes", "features", "on_disc_matrix", "on_disc_matrix")), if (is.null(idx)) "" else paste0("_", idx), c(".mtx", ".tsv", ".tsv", ".rds", ".h5"))
  to_save_locs <- setNames(paste0(data_dir, "/", f_names),  c("mtx", "barcodes", "features", "on_disc_matrix",  "on_disc_matrix_h5"))
  return(to_save_locs)
}


#' Get metadata odm list
#'
#' @param mat_list a list of matrices
#' @param idx_start the starting index
#' @param temp_test_dir directory used by the tests
#'
#' @return a list of initialized metadata_odms
#' @noRd
get_metadata_odm_list <- function(mat_list, idx_start, temp_test_dir) {
  cov_odms <- vector(mode = "list", length = length(mat_list))
  for (i in seq(1, length(mat_list))) {
    fps <- get_simulation_data_fps(data_dir = temp_test_dir, idx = i + idx_start - 1L)
    # check if on_disc_matrix already has been created; if so, delete that as well as the .h5 file
    if (file.exists(fps[["on_disc_matrix_h5"]])) file.remove(fps[["on_disc_matrix_h5"]]) %>% invisible() # h5 file
    m <- mat_list[[i]]
    n_data_points <- length(m@x)
    if (stats::rbinom(1, 1, 0.5)) {
      chunk_size <- sample(x = seq(2, n_data_points - 1), size = 1) # choose chunk size less than n_data_points
    } else {
      chunk_size <- sample(x = seq(n_data_points + 1, 2 * n_data_points), size = 1)
    }
    cov_odm_obj <- create_ondisc_matrix_from_mtx(mtx_fp = fps[["mtx"]],
                                                 barcodes_fp = fps[["barcodes"]],
                                                 features_fp = fps[["features"]],
                                                 n_lines_per_chunk = chunk_size,
                                                 on_disk_dir = temp_test_dir,
                                                 return_metadata_ondisc_matrix = TRUE)
    cov_odms[[i]] <- cov_odm_obj
  }
  return(cov_odms)
}
