#' Return sparse matrix from contiguous index
#'
#' @param h5_file the file path to the h5 file
#' @param contiguous_idx_range a length-two INTEGER vector containing the starting end ending index (inclusive)
#' @param index_on_cell (boolean) index on cell (TRUE) or feature (FALSE)
#' @param logical_mat (boolean)
#'
#' @return a Matrix object with the given contiguous_idx_range.
return_spMatrix_contig_index <- function(h5_file, contiguous_idx_range, index_on_cell, logical_mat) {
  p_name <- paste0(if(index_on_cell) "cell" else "feature" , "_ptr")
  idx_name <- paste0(if(index_on_cell) "feature" else "cell", "_idxs")
  if (!logical_mat) dat_name <- paste0("data_", if(index_on_cell) "csc" else "csr")
  start <- contiguous_idx_range[1]
  count <- contiguous_idx_range[2] - start + 2L # add 2 regardless of 0/1-based indexing
  # extract pointer and center at zero.
  p <- as.integer(rhdf5::h5read(file = h5_file, name = p_name, start = start, count = count) + 1L)
  p_start <- p[1]
  centered_p <- p[-1] - p[1]
  p_count <- centered_p[count - 1L]
  if (p_count > 0) { # if nonempty
    idxs <- as.integer(rhdf5::h5read(file = h5_file, name = idx_name, start = p_start, count = p_count))
    if (!logical_mat) {
      dat <- as.numeric(rhdf5::h5read(file = h5_file, name = dat_name, start = p_start, count = p_count))
    }
  } else { # if empty
    idxs <- integer(0)
    if (!logical_mat) dat <- numeric(0)
  }
  # determine dimension of submatrix; initialize full pointer
  full_p <- c(0L, centered_p)
  dimen <- rhdf5::h5read(file = h5_file, name = "dimension") %>% as.integer()
  dim_to_pass <- if (index_on_cell) c(dimen[1], count - 1L) else c(count - 1L, dimen[2])
  # finally, create matrix
  ret <- if (index_on_cell && !logical_mat) { # index on cell, integer matrix
    new(getClass(Class = "dgCMatrix", where = "Matrix"),
        x = dat,
        i = idxs,
        p = full_p,
        Dim = dim_to_pass)
  } else if (!index_on_cell && !logical_mat) { # index on feature, integer matrix
    new(getClass(Class = "dgRMatrix", where = "Matrix"),
        x = dat,
        j = idxs,
        p = full_p,
        Dim = dim_to_pass)
  } else if (index_on_cell && logical_mat) { # index on cell, logical matrix
    new(getClass(Class = "lgCMatrix", where = "Matrix"),
        x = rep(TRUE, length(idxs)),
        i = idxs,
        p = full_p,
        Dim = dim_to_pass)
  } else { # index on feature, logical matrix
    new(getClass(Class = "lgRMatrix", where = "Matrix"),
        x = rep(TRUE, length(idxs)),
        j = idxs,
        p = full_p,
        Dim = dim_to_pass)
  }
  return(ret)
}


#' return_spMatrix_from_index
#'
#' @param h5_file an h5 file
#' @param subset_vector an ordered, possibility discontiguous set of indexes
#' @param index_on_cell (boolean) index on cell (TRUE) or feature (FALSE)
#' @param logical_mat (boolean) is matrix logical (FALSE) or integer (TRUE)
#'
#' @return a sparse matrix subset according to idx
return_spMatrix_from_index <- function(h5_file, subset_vector, index_on_cell, logical_mat) {
  if (length(subset_vector) == 1) {
    ret <- return_spMatrix_contig_index(h5_file, c(subset_vector, subset_vector), index_on_cell, logical_mat)
  }
  intervals <- find_contig_subseqs(subset_vector)
  start_ints <- intervals[[1]]
  end_ints <- intervals[[2]]
  n_intervals <- length(start_ints)
  if (n_intervals == 1) {
    ret <- return_spMatrix_contig_index(h5_file, c(start_ints, end_ints), index_on_cell, logical_mat)
  } else {
    mats <- lapply(X = seq(1, n_intervals), FUN = function(i) {
      curr_int <- c(start_ints[i], end_ints[i])
      return_spMatrix_contig_index(h5_file, curr_int, index_on_cell, logical_mat)
    })
    ret <- do.call(if(index_on_cell) cbind else rbind, mats)
  }
  return(ret)
}
