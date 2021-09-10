# nocov start
# These are helper functions used for testing purposes only. None of these functions are exported.

#' Create new directory
#'
#' Creates a unique directory; returns a file path to that directory.
#'
#' @return file path to a unique, empty directory
#' @export
create_new_directory <- function() {
  f <- tempfile()
  while (dir.exists(f)) f <- tempfile()
  dir.create(f)
  return(f)
}


#' Create a random matrix
#'
#' Creates and returns a random matrix.
#'
#' @param n_row number of rows
#' @param n_col number of columns
#' @param logical_mat boolean indicating whether the matrix is logical (TRUE) or integer (FALSE)
#' @param p_zero probability an entry will be zero
#' @param p_set_col_zero fraction of columns to set to zero
#' @param p_set_row_zero fraction of rows to set to zero
#' @param matrix_values set of values from which to draw the matrix entries (applicable to integer matrices only)
#'
#' @return a randomly-generated matrix in sparse format
create_random_matrix <- function(n_row, n_col, logical_mat, p_zero = 0.95, p_set_col_zero = 0.05, p_set_row_zero = 0.05, matrix_values = seq(1L, 10L)) {
  # sample binary matrix
  r <- matrix(data = stats::rbinom(n =  n_row * n_col, size = 1, prob = 1 - p_zero), nrow = n_row, ncol = n_col)
  # randomly set some rows or columns to all zero
  zero_rows <- sample(seq(1, n_row), size = floor(p_set_row_zero * n_row), replace = FALSE)
  zero_cols <- sample(seq(1, n_col), size = floor(p_set_row_zero * n_col), replace = FALSE)
  r[zero_rows,] <- 0; r[,zero_cols] <- 0
  # initialize integer or logical matrix
  if (!logical_mat) {
    m <- matrix(data = sample(x = matrix_values, size = n_row * n_col, replace = TRUE), nrow = n_row, ncol = n_col)
    out <- m * r
  } else {
    out <- r == 1
  }
  # convert the matrix to sparse format and return
  out <- Matrix::Matrix(data = out, sparse = TRUE)
  return(out)
}


#' Save random matrix as a 10X object
#'
#' Writes an R matrix to disk in .mtx format. Writes barcodes.tsv and features.tsv files as well. genes.
#'
#' @param m a sparse Matrix object
#' @param data_dir the directory in which to store the matrix
#' @param cell_barcodes the cell barcodes
#' @param feature_df data frame of features
#'
#' @return the file paths to the matrix.mtx, barcodes.tsv, and features.tsv files.
#' @noRd
save_random_matrix_as_10x <- function(m, data_dir, cell_barcodes, feature_df) {
  fps <- paste0(data_dir, "/", c("matrix.mtx", "barcodes.tsv", "features.tsv"))
  names(fps) <- c("matrix", "barcodes", "features")
  Matrix::writeMM(obj = m, file = fps["matrix"])
  readr::write_tsv(x = data.frame(cell_barcodes), file = fps["barcodes"], col_names = FALSE)
  readr::write_tsv(x = feature_df, file = fps["features"], col_names = FALSE)
  return(fps)
}


#' Save random matrix as a .h5 file
#'
#' Writes an R matrix to disk in .h5 format.
#'
#' @param m a sparse Matrix object
#' @param data_dir the directory in which to store the matrix
#' @param cell_barcodes the cell barcodes
#' @param feature_df data frame of features
#'
#' @return the file paths to the matrix.mtx, barcodes.tsv, and features.tsv files.
#' @noRd
save_random_matrix_as_h5 <- function(m, data_dir, cell_barcodes, feature_df) {
  h5_fp <- paste0(data_dir, "/matrix.h5")
  csc_r_matrix <- as(m, "dgCMatrix")
  # Write CSC
  rhdf5::h5write(csc_r_matrix@p, file = h5_fp, name = "indptr")
  rhdf5::h5write(csc_r_matrix@i, file = h5_fp, name = "indices")
  rhdf5::h5write(csc_r_matrix@x, file = h5_fp, name = "data")
  rhdf5::h5write(csc_r_matrix@Dim, file = h5_fp, name = "shape")
  rhdf5::h5write(cell_barcodes, file = h5_fp, name = "barcodes")
  rhdf5::h5write(feature_df$gene_ids, file = h5_fp, name = "genes")
  rhdf5::h5write(feature_df$gene_names, file = h5_fp, name = "gene_names")
  return(h5_fp)
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
#' @noRd
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
#' Creates synthetic data for testing purposes. Optionally writes the data to disk.
#'
#' @param n_row number of rows
#' @param n_col number of columns
#' @param logical_mat boolean indicating whether a given dataset is logical (TRUE) or integer (FALSE).
#' @param start_pos position to start numbering the barcodes
#' @param write_as_mtx_to_disk boolean indicating whether to write the mtx file, barcodes file, and features file to disk as side-effect.
#' @param write_as_h5_to_disk boolean indicating whether to write the .h5 list of files to disk as side-effect.
#' @param file_dir directory in which to write the data to disk
#' @return a list:
#' (i) R matrix
#' (ii) in-memory cell barcodes
#' (iii) in-memory features df
#' (iv) number of nonzero entries in expression matrix
#' (v) file path to the .mtx file
#' (vi) file path to the barcodes.tsv file
#' (vii) file path to the features.tsv file
#' (viii) OR a file path to the .h5 file
#' @export
#'
#' @examples
#' n_row <- 1000
#' n_col <- 500
#' logical_mat <- FALSE
#' synth_data <- create_synthetic_data(n_row, n_col, logical_mat)
#' # side-effect: creates .mtx file, barcodes.tsv file, and features.tsv file on disk.
create_synthetic_data <- function(n_row, n_col, logical_mat, start_pos = 0L, write_as_mtx_to_disk = TRUE, write_as_h5_to_disk = FALSE, file_dir = NULL) {
  # first, create an in-memory r matrix
  in_mem_r_mat <- create_random_matrix(n_row, n_col, logical_mat)
  # next, create a character vector of barcodes
  barcodes <- paste0("cell_", seq(1+start_pos, n_col+start_pos))
  # next, create features dfs
  gene_names <- paste0("gene_", seq(1, n_row))
  idxs <- sort(sample(x = seq(1, n_row), size = floor(n_row/10), replace = FALSE))
  gene_names[idxs] <- paste0("MT-", idxs)
  gene_ids <- paste0("ENSG000", seq(1, n_row))
  features_df <- data.frame(gene_ids, gene_names)
  # create "out," containing the matrix, gene ids, and features df
  out <- list(r_mat = in_mem_r_mat, barcodes = barcodes,
              features_df = features_df, n_nonzero = length(in_mem_r_mat@x))
  if (write_as_h5_to_disk) {
    # create the directory
    if (is.null(file_dir)) file_dir <- create_new_directory()
    h5_fp <- save_random_matrix_as_h5(m = in_mem_r_mat,
                                     data_dir = file_dir,
                                     cell_barcodes = barcodes,
                                     feature_df = features_df)
    out$h5_fp <- h5_fp
  }
  if (write_as_mtx_to_disk) {
    # create the directory
    if (is.null(file_dir)) file_dir <- create_new_directory()
    fps <- save_random_matrix_as_10x(m = in_mem_r_mat,
                                     data_dir = file_dir,
                                     cell_barcodes = barcodes,
                                     feature_df = features_df)
    out$matrix_fp <- fps[["matrix"]]
    out$barcodes_fp <- fps[["barcodes"]]
    out$features_fp <- fps[["features"]]
  }
  return(out)
}


#' get_random_subset
#'
#' Return a random integer subset (of random length) of the set [n].
#'
#' @param n an integer
#'
#' @return a random subset of 1, 2, ..., n.
get_random_subset <- function(n) {
  sample(x = seq(1, n), size = sample(x = seq(1, n), size = 1L), replace = FALSE)
}

# nocov end
