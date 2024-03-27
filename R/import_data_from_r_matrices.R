#' Create ODM from R matrix
#'
#' `create_odm_from_r_matrix()` creates an `odm` object from an input R matrix.
#'
#' @note
#' Users do not need to specify the `integer_id` argument; this argument is used by functions within the `sceptre` package that call `create_odm_from_r_matrix()`.
#'
#' @param mat an R matrix of class "matrix", "dgCMatrix", "dgRMatrix", or "dgTMatrix"
#' @param file_to_write file path specifying the location in which to write the backing `.odm` file
#' @param chunk_size chunk size to use in the backing HDF5 file
#' @param compression_level compression level to use in the backing HDF5 file
#' @param integer_id integer ID to write to the backing `.odm` file; users need not specify this argument.
#'
#' @return an odm object representing the data
#'
#' @examples
#' library(Matrix)
#' set.seed(4)
#' n_row <- 500
#' n_col <- 800
#' v <- rpois(n = n_row * n_col, lambda = 1)
#' v[rbinom(n = n_row * n_col, size = 1, prob = 0.8) == 1] <- 0L
#' mat <- matrix(data = v, nrow = n_row, ncol = n_col)
#' rownames(mat) <- paste0("gene_", seq_len(nrow(mat)))
#' file_to_write <- paste0(tempdir(), "/gene.odm")
#' odm <- create_odm_from_r_matrix(mat, file_to_write)
create_odm_from_r_matrix <- function(mat, file_to_write, chunk_size = 1000L, compression_level = 3L, integer_id = 0L) {
  # convert the matrix into a dgRMatrix using the sceptre function
  mat <- sceptre:::set_matrix_accessibility(mat)

  # expand tilde for the file to write
  file_to_write <- expand_tilde(file_to_write)

  # check for rownames
  if (is.null(rownames(mat))) {
    stop("Row names must be supplied for `mat`.")
  }

  # verify that chunk_size does not exceed data size
  if (length(mat@j) <= chunk_size) {
    stop("Decrease the chunk size.")
  }

  # initialize odm
  create_odm_r_matrix_cpp(file_name_in = file_to_write,
                          feature_ids = rownames(mat),
                          n_features = nrow(mat),
                          n_cells = ncol(mat),
                          integer_id = integer_id,
                          chunk_size = chunk_size,
                          compression_level = compression_level,
                          j = mat@j, x = as.integer(mat@x), p = mat@p)

  # create odm object
  out <- initialize_odm_from_backing_file(odm_file = file_to_write)
  return(out)
}
