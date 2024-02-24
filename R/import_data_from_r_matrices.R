#' Create ODM from R matrix
#'
#' @param mat an R matrix of class "dgRMatrix"
#' @param file_to_write name of the file in which to write the backing `.odm` file
#' @param chunk_size chunk size to use in the backing HDF5 file
#' @param compression_level compression level to use in the backing HDF5 file
#' @param integer_id ID used to identify the `.odm` file
#'
#' @return an odm object representing the data
#' @export
#'
#' @examples
#' set.seed(4)
#' n_row <- 500
#' n_col <- 800
#' v <- rpois(n = n_row * n_col, lambda = 1)
#' v[rbinom(n = n_row * n_col, size = 1, prob = 0.8) == 1] <- 0L
#' m <- matrix(data = v, nrow = n_row, ncol = n_col)
#' rownames(m) <- paste0("gene_", seq(1L, n_row))
#' mat <- as(m, "RsparseMatrix")
#' odm <- create_odm_from_r_matrix(mat, paste0(tempdir(), "gene.odm"))
create_odm_from_r_matrix <- function(mat, file_to_write, chunk_size = 1000L, compression_level = 3L, integer_id = 0L) {
  # 0. check that the class is correct
  if (!is(mat, "dgRMatrix")) stop("`mat` must be an object of class `dgRMatrix`.")
  file_to_write <- expand_tilde(file_to_write)
  if (is.null(directory_to_write)) {
    stop("`directory_to_write` cannot be `NULL`.")
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
