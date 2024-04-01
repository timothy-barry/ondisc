#' Create ODM from R matrix
#'
#' `create_odm_from_r_matrix()` converts an R matrix (stored in standard dense format or sparse format) into an `odm` object.
#'
#' @param mat an R matrix of class `"matrix"`, `"dgCMatrix"`, `"dgRMatrix"`, or `"dgTMatrix"`.
#' @param file_to_write a fully-qualified file path specifying the location in which to write the backing `.odm` file.
#' @param chunk_size (optional; default `1000L`) a positive integer specifying the chunk size to use to store the data in the backing HDF5 file.
#' @param compression_level (optional; default `3L`) an integer in the inveral \[0, 9\] specifying the compression level to use to store the data in the backing HDF5 file.
#'
#' @return an `odm` object
#' @export
#'
#' @examples
#' library(sceptredata)
#' data(lowmoi_example_data)
#' gene_matrix <- lowmoi_example_data$response_matrix
#' file_to_write <- paste0(tempdir(), "gene.odm")
#' odm_object <- create_odm_from_r_matrix(
#'   mat = gene_matrix,
#'   file_to_write = file_to_write
#' )
create_odm_from_r_matrix <- function(mat, file_to_write, chunk_size = 1000L, compression_level = 3L) {
  create_odm_from_r_matrix_internal(mat = mat, file_to_write = file_to_write,
                                    chunk_size = chunk_size, compression_level = compression_level)
}


create_odm_from_r_matrix_internal <- function(mat, file_to_write, chunk_size = 1000L, compression_level = 3L, integer_id = 0L) {
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
