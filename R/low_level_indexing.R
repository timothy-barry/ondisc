#' Return sparse matrix from contiguous index
#'
#' Takes an index range (length two integer vector), h5_file, and h5_group (column or row) and returns the associated data, indexes, and pointer.
#'
#' @param h5_file path to h5 file on disk.
#' @param h5_group h5 group from which to extract data (compressed_sparse_row/ or compressed_sparse_column/)
#' @param contiguous_idx_range an integer vector of length 2 specifying the (inclusive) range of rows/columns to extract
#'
#' @return a list containing the 0-based pointer vector (excluding 0 itself) and the associated data
return_spMatrix_from_contiguous_index <- function(h5_file, h5_group, contiguous_idx_range) {
  p_range <- rhdf5::h5read(file = h5_file, name = paste0(h5_group, "p"), index = list(contiguous_idx_range[1]:(contiguous_idx_range[length(contiguous_idx_range)] + 1)))
  n_vectors_p1 <- length(p_range)
  if (p_range[n_vectors_p1] - 1 >= p_range[1]) {
    data <- rhdf5::h5read(file = h5_file, name = paste0(h5_group, "data"), index = list(p_range[1]:(p_range[n_vectors_p1] - 1), NULL))
    p_out <- p_range[-1] - p_range[1]
  } else {
    data <- matrix(data = integer(), ncol = 2, nrow = 0)
    p_out <- rep(0, n_vectors_p1 - 1)
  }
  list(p = p_out, data = data)
}


#' Return sparse matrix from index
#'
#' Given an h5 file containing an on_disc_matrix, returns the sparse matrix corresponding to an index subset.
#'
#' @param h5_file path to h5 file on disk
#' @param idx indexes to extract
#' @param col_idx boolean; extract column indexes (true) or row indexes (false)
#'
#' @export
#' @return a sparse matrix of class Matrix.
#' @examples
#' h5_file <- system.file("extdata", "on_disc_matrix_1.h5", package = "ondisc")
#' if (h5_file != "") {
#' idx <- 2:3
#' col_idx <- TRUE
#' return_spMatrix_from_index(h5_file, idx, col_idx)
#' }
return_spMatrix_from_index <- function(h5_file, idx, col_idx) {
  # assume idx is a list of non-negative, unique integers. First, determine if idx is sorted.
  idx_unsorted <- is.unsorted(idx)
  if (idx_unsorted) permutation <- order(x = idx, method = "radix")

  # Next, determine the distinct subsequences of contiguous integers
  contiguous_intervals <- R.utils::seqToIntervals(idx)
  h5_group <- paste0("compressed_sparse_", if (col_idx) "column" else "row",  "/")
  n_contiguous_intervals <- nrow(contiguous_intervals)

  # If there are multiple contiguous intervals, then load each separately and combine.
  if (n_contiguous_intervals >= 2) {
    sub_matrices <- purrr::map(.x = 1:n_contiguous_intervals, .f = function(i) return_spMatrix_from_contiguous_index(h5_file, h5_group, contiguous_intervals[i,]))
    full_matrix_data <- purrr::map(sub_matrices, function(i) i[["data"]]) %>% do.call(what = rbind, args = .)
    sub_pointers <- purrr::map(sub_matrices, function(i) i[["p"]])
    for (i in 2:(n_contiguous_intervals)) {
      sub_pointers[[i]] <- sub_pointers[[i]] + sub_pointers[[i - 1]][length(sub_pointers[[i - 1]])]
    }
    full_p <- c(0, unlist(sub_pointers)) %>% as.integer()
  } else { # If there is only a single contiguous interval, then load that only.
    mat <- return_spMatrix_from_contiguous_index(h5_file, h5_group, contiguous_intervals)
    full_matrix_data <- mat$data
    full_p <- c(0, mat$p) %>% as.integer()
  }
  Dim <- rhdf5::h5read(file = h5_file, name = "metadata/dimension") %>% as.integer()
  dim_to_pass <- if (col_idx) c(Dim[1], length(idx)) else c(length(idx), Dim[2])
  if (col_idx) {
    ret <- new(Class = "dgCMatrix", x = as.numeric(full_matrix_data[,2]), i = full_matrix_data[,1], p = full_p, Dim = dim_to_pass)
  } else {
    ret <- new(Class = "dgRMatrix", x = as.numeric(full_matrix_data[,2]), j = full_matrix_data[,1], p = full_p, Dim = dim_to_pass)
  }
  # Finally, if the indexes were unsorted, reorder the matrix.
  if (idx_unsorted) ret <- if (col_idx) ret[,Matrix::invPerm(permutation)] else ret[Matrix::invPerm(permutation),]
  return(ret)
}
