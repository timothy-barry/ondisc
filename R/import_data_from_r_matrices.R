#' Create `odm` object from R matrix
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
#' file_to_write <- paste0(tempdir(), "/gene.odm")
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

# helper functions, copied and pasted from sceptre
set_matrix_accessibility <- function(matrix_in, make_row_accessible = TRUE) {
  # if a logical matrix, convert to the corresponding numeric matrix; consider more efficient implementation later
  if (methods::is(matrix_in, "lgRMatrix")) {
    attr(matrix_in, "class") <- "dgRMatrix"
    matrix_in@x <- rep(1.0, length(matrix_in@j))
  }
  if (methods::is(matrix_in, "lgCMatrix")) {
    attr(matrix_in, "class") <- "dgCMatrix"
    matrix_in@x <- rep(1.0, length(matrix_in@i))
  }
  if (methods::is(matrix_in, "lgTMatrix")) {
    attr(matrix_in, "class") <- "dgTMatrix"
    matrix_in@x <- rep(1.0, length(matrix_in@j))
  }

  if (methods::is(matrix_in, "dgRMatrix") && make_row_accessible) {
    out <- matrix_in
  } else if (methods::is(matrix_in, "dgCMatrix") && !make_row_accessible) {
    out <- matrix_in
  } else if (methods::is(matrix_in, "odm")) {
    out <- matrix_in
  } else {
    # basic info: get matrix_in, n_responses, response_ids, and n_cells
    n_responses <- nrow(matrix_in)
    response_ids <- rownames(matrix_in)
    n_cells <- ncol(matrix_in)
    if (methods::is(matrix_in, "dgTMatrix")) {
      i <- matrix_in@i
      j <- matrix_in@j
    } else if (methods::is(matrix_in, "matrix")) {
      matrix_in <- methods::as(matrix_in, "TsparseMatrix")
      i <- matrix_in@i
      j <- matrix_in@j
    } else if (methods::is(matrix_in, "dgCMatrix")) {
      i <- matrix_in@i
      j <- convert_pointer_to_index_vector_v2(matrix_in@p)
    } else if (methods::is(matrix_in, "dgRMatrix")) {
      i <- convert_pointer_to_index_vector_v2(matrix_in@p)
      j <- matrix_in@j
    } else {
      stop("Class not recognized.")
    }
    x <- matrix_in@x

    # perform the sorting operation
    sort_in_place <- methods::is(matrix_in, "matrix") || methods::is(matrix_in, "dgTMatrix")
    if (sort_in_place) {
      l <- list(i = i, j = j, x = x)
      dt <- data.table::setDT(l) # assign by reference
    } else {
      dt <- data.table::data.table(i = i, j = j, x = x) # deep copy
    }

    # finally, initialize the output
    if (make_row_accessible) {
      data.table::setorderv(x = dt, cols = c("i", "j"))
      p <- obtain_pointer_vector(i = dt$i, dim = n_responses)
      out <- Matrix::sparseMatrix(j = 1, p = c(0, 1), x = 1, repr = "R")
      out@j <- dt$j
    } else {
      data.table::setorderv(x = dt, cols = c("j", "i"))
      p <- obtain_pointer_vector(i = dt$j, dim = n_cells)
      out <- Matrix::sparseMatrix(i = 1, p = c(0, 1), x = 1, repr = "C")
      out@i <- dt$i
    }
    out@p <- p
    out@x <- dt$x
    out@Dim <- c(n_responses, n_cells)
    rownames(out) <- response_ids
  }
  return(out)
}


convert_pointer_to_index_vector_v2 <- function(p) {
  l <- diff(p)
  rep(x = seq(0, length(l) - 1L), times = l)
}
