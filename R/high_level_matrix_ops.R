#' Get on_disc subset vector
#'
#' Returns the subset vector (either feature or cell) for an ondisc_matrix object.
#'
#' @param x the ondisc_matrix.
#' @param cell_idx (boolean) TRUE = cell, FALSE = feature.
#'
#' @return the requested subset indexes
#' @noRd
get_subset_vector <- function(x, cell_idx) {
  idx_slot <- paste0(if (cell_idx) "cell" else "feature", "_subset")
  return(slot(x, idx_slot))
}


#' Get dim
#'
#' Returns the dimension of an on_disc matrix.
#'
#' @param x an ondisc_matrix object
#'
#' @return the dimension
#' @noRd
get_dim <- function(x) {
  feature_unsubset <- identical(x@feature_subset, NA_integer_)
  cell_unsubset <- identical(x@cell_subset, NA_integer_)
  underlying_dim <- x@underlying_dimension
  out <- c(if (feature_unsubset) underlying_dim[1] else length(x@feature_subset),
           if (cell_unsubset) underlying_dim[2] else length(x@cell_subset))
  return(out)
}


#' Get names helper function
#'
#' @param x an ondisc_matrix
#' @param name_to_get the names to extract from x (one of cell_barcodes, feature_ids, feature_names)
#'
#' @return a character vector containing the requested names
#' @noRd
get_names <- function(x, name_to_get) {
  return(slot(x, name_to_get))
}


#' Subset by feature or cell
#'
#' Subsets an ondisc_matrix by either feature or cell.
#'
#' @param x an ondisc_matrix
#' @param idx a numeric, character, or logical index
#' @param subset_on_cell (boolean) subset on cell (TRUE) or feature (FALSE)
#' @param subset_string_arrays (boolean) subset the character vectors (i.e., cell barcodes, feature IDs, and feature names) as well? This is necessary only when returning a subset matrix, not when preparing to extract a submatrix into memory.
#'
#' @return a subset ondisc_matrix
#' @export
subset_by_feature_or_cell <- function(x, idx, subset_on_cell, subset_string_arrays = FALSE) {
  subset_slot <- paste0(if (subset_on_cell) "cell" else "feature", "_subset")
  if (identical(slot(x, subset_slot), NA_integer_)) {
    n_elements <- if (subset_on_cell) ncol(x) else nrow(x)
    slot(x, subset_slot) <- seq(1, n_elements)
  }
  # perform different subset operation based on class of idx
  if (is(idx, "numeric")) {
    if (max(idx) > length(slot(x, subset_slot))) stop("Numeric index out of bounds.")
  } else if (is(idx, "logical")) {
    if (length(idx) > length(slot(x, subset_slot))) stop("Logical vector too long.")
  } else if (is(idx, "character")) {
    curr_ids <- slot(x, if (subset_on_cell) "cell_barcodes" else "feature_ids")
    # update to integer idx
    idx <- match(x = idx, table = curr_ids)
    if (any(is.na(idx))) stop(paste(if (subset_on_cell) "Cell barcode" else "feature id","not present in data."))
  } else {
    stop("idx must be of type numeric, character, or logical.")
  }
  # perform the subset on the subset_slot
  slot(x, subset_slot) <- slot(x, subset_slot)[idx]
  # Also index into the string arrays, if necessary.
  if (subset_string_arrays) {
    if (subset_on_cell) {
      slot(x, "cell_barcodes") <- slot(x, "cell_barcodes")[idx]
    } else {
      slot(x, "feature_names") <- slot(x, "feature_names")[idx]
      slot(x, "feature_ids") <- slot(x, "feature_ids")[idx]
    }
  }
  return(x)
}


#' Extract matrix
#'
#' @param x an ondisc_matrix.
#'
#' @return an in-memory version of x in the form of a Matrix object.
#' @export
extract_matrix <- function(x) {
  # First, determine which axis to index on; always index on the shorter axis.
  x_dim <- dim(x)
  if (x_dim[1] == 0 || x_dim[2] == 0) {
    stop("ondisc_matrix has 0 rows or columns; cannot extract sub-matrix into memory.")
  }
  index_on_cell <- x_dim[2] < x_dim[1]
  subset_vector <- get_subset_vector(x, index_on_cell)
  # if subset_vector is NA, initialize it as an integer vector (this will happen only rarely)
  if (identical(subset_vector, NA_integer_)) {
    subset_vector <- seq(1, if (index_on_cell) x_dim[2] else x_dim[1])
  }
  # return submatrix
  out <- return_spMatrix_from_index(x@h5_file, subset_vector, index_on_cell, x@logical_mat, x@underlying_dimension)
  # If necessary, subset along the other axis
  second_subset <- get_subset_vector(x, !index_on_cell)
  if (!identical(second_subset, NA_integer_)) {
    out <- if (index_on_cell) out[second_subset,,drop=FALSE] else out[,second_subset,drop=FALSE]
  }
  return(out)
}
