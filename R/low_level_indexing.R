#' Return sparse matrix from index
#'
#' @param h5_file the file path to the h5 file
#' @param subset_vector an integer vector giving the subset
#' @param index_on_cell (boolean) index on cell (TRUE) or feature (FALSE)
#' @param logical_mat (boolean)
#' @param underlying_dimension a length 2 integer vector giving the underlying dimension of the expression matrix.
#'
#' @return a Matrix object.
return_spMatrix_from_index <- function(h5_file, subset_vector, index_on_cell, logical_mat, underlying_dimension) {
  # define variables specific to index_on_cell and logical_mat
  p_name <- paste0(if (index_on_cell) "cell" else "feature" , "_ptr")
  idx_name <- paste0(if (index_on_cell) "feature" else "cell", "_idxs")
  umi_counts_name <- paste0("data_", if (index_on_cell) "csc" else "csr")
  dim_to_pass <- if (index_on_cell) {
    c(underlying_dimension[1], length(subset_vector))
  } else {
    c(length(subset_vector), underlying_dimension[2])
  }
  subset_vector <- subset_vector - 1
  # get the matrix
  flat_mat <- index_h5_file(file_name_in = h5_file, p_name_in = p_name, idx_name_in = idx_name, umi_counts_name_in = umi_counts_name, subset_vector = subset_vector)
  p <- flat_mat[[1]]
  idxs <- flat_mat[[2]]
  dat <- as.numeric(flat_mat[[3]])

  # finally, create the matrix
  ret <- if (index_on_cell && !logical_mat) { # index on cell, integer matrix
    new(getClass(Class = "dgCMatrix", where = "Matrix"),
        x = dat,
        i = idxs,
        p = p,
        Dim = dim_to_pass)
  } else if (!index_on_cell && !logical_mat) { # index on feature, integer matrix
    new(getClass(Class = "dgRMatrix", where = "Matrix"),
        x = dat,
        j = idxs,
        p = p,
        Dim = dim_to_pass)
  } else if (index_on_cell && logical_mat) { # index on cell, logical matrix
    new(getClass(Class = "lgCMatrix", where = "Matrix"),
        x = rep(TRUE, length(idxs)),
        i = idxs,
        p = p,
        Dim = dim_to_pass)
  } else { # index on feature, logical matrix
    new(getClass(Class = "lgRMatrix", where = "Matrix"),
        x = rep(TRUE, length(idxs)),
        j = idxs,
        p = p,
        Dim = dim_to_pass)
  }
  return(ret)
}
