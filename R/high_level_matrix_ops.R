#' Get on_disc subset vector
#'
#' Returns the subset vector and order vector (either feature or cell) for an ondisc_matrix object.
#'
#' @param x the ondisc_matrix.
#' @param cell_idx (boolean) TRUE = cell, FALSE = feature.
#'
#' @return the requested subset indexes
get_subset_vector <- function(x, cell_idx) {
  idx_slot <- paste0(if (cell_idx) "cell" else "feature", "_subset")
  return(slot(x, idx_slot))
}

#' @rdname get_subset_vector
get_order_vector <- function(x, cell_idx) {
  idx_slot <- paste0(if (cell_idx) "cell" else "feature", "_subset_order")
  return(slot(x, idx_slot))
}

#' Get dim
#'
#' Returns the dimension of an on_disc matrix.
#'
#' @param x an ondisc_matrix object
#'
#' @return the dimension
get_dim <- function(x) {
  cell_unsubset <- identical(x@cell_subset, NA_integer_)
  feature_unsubset <- identical(x@feature_subset, NA_integer_)

  if (cell_unsubset || feature_unsubset) dim_h5 <- rhdf5::h5read(file = x@h5_file, name = "dimension") %>% as.integer()
  n_row_out <- if (feature_unsubset) dim_h5[1] else length(x@feature_subset)
  n_col_out <- if (cell_unsubset) dim_h5[2] else length(x@cell_subset)
  return(c(n_row_out, n_col_out))
}


#' Get names helper function
#'
#' @param x an ondisc_matrix
#' @param name_to_get the names to extract from x (one of cell_barcodes, feature_ids, feature_names)
#'
#' @return a character vector containing the requested names
get_names <- function(x, name_to_get) {
  names_out <- rhdf5::h5read(file = x@h5_file, name = name_to_get) %>% as.character()
  idx <- get_subset_vector(x, name_to_get == "cell_barcodes")
  if (identical(idx, NA_integer_)) names_out else names_out[idx]
}


#' Get feature names, feature IDs, and cell barcodes
#'
#' Functions to extract cell barcodes, feature ids, and feature names from an ondisc_matrix.
#'
#' @param x An ondisc_matrix.
#' @return The cell barcodes of this on-disc matrix.
#' @export
#' @examples
#' # NOTE: You must create the HDF5 file "example.h5" to run this example.
#' # Navigate to the help file of "create_ondisc_matrix_from_10x_mtx"
#' # (via ?create_ondisc_matrix_from_10x_mtx), and execute the code in the example.
#' odm_fp <- system.file("extdata", "example.h5", package = "ondisc")
#' if (odm_fp != "") { # if required file exists, ...
#' odm <- ondisc_matrix(h5_file = odm_fp)
#' barcodes <- get_cell_barcodes(odm)
#' feature_ids <- get_feature_ids(odm)
#' feature_names <- get_feature_names(odm)
#' }
get_cell_barcodes <- function(x) {
  get_names(x, "cell_barcodes")
}


#' @rdname get_cell_barcodes
#' @param x An ondisc_matrix.
#' @export
#' @return The cell feature ids of this on-disc matrix.
get_feature_ids <- function(x) {
  get_names(x, "feature_ids")
}


#' @rdname get_cell_barcodes
#' @param x An ondisc_matrix.
#' @export
#' @return The feature names of this on-disc matrix.
get_feature_names <- function(x) {
  get_names(x, "feature_names")
}


#' Subset by feature or cell
#'
#' Subsets an ondisc_matrix by either feature or cell.
#'
#' @param x an ondisc_matrix
#' @param idx a numeric, character, or logical index
#' @param subset_on_cell (boolean) subset on cell (TRUE) or feature (FALSE)
#'
#' @return a subset ondisc_matrix
subset_by_feature_or_cell <- function(x, idx, subset_on_cell) {
  subset_slot <- paste0(if (subset_on_cell) "cell" else "feature", "_subset")
  order_slot <- paste0(subset_slot, "_order")
  if (identical(slot(x, subset_slot), NA_integer_)) {
    n_elements <- if (subset_on_cell) ncol(x) else nrow(x)
    slot(x, subset_slot) <- seq(1, n_elements)
  }
  # perform different subset operation based on class of idx
  if (is(idx, "numeric")) {
    if (max(idx) > length(slot(x, subset_slot))) stop("Numeric index out of bounds.")
    slot(x, subset_slot) <- slot(x, subset_slot)[idx]
  } else if (is(idx, "logical")) {
    if (length(idx) > length(slot(x, subset_slot))) stop("Logical vector too long.")
    slot(x, subset_slot) <- slot(x, subset_slot)[idx]
  } else if (is(idx, "character")) {
    all_ids <- rhdf5::h5read(file = x@h5_file, name = if (subset_on_cell) "cell_barcodes" else "feature_ids")
    curr_ids <- all_ids[slot(x, subset_slot)]
    within_cur_ids_idxs <- match(x = idx, table = curr_ids)
    if (any(is.na(within_cur_ids_idxs))) stop(paste(if(subset_on_cell) "Cell barcode" else "feature id","not present in data."))
    slot(x, subset_slot) <- slot(x, subset_slot)[within_cur_ids_idxs]
  } else {
    stop("idx must be of type numeric, character, or logical.")
  }
  # check for duplicates and store ordering
  l <- length(slot(x, subset_slot))
  ordering <- order(slot(x, subset_slot), method = "radix")
  tmp <- rle(slot(x, subset_slot)[ordering])
  if (any(tmp$lengths >= 2L)) stop("Duplicate indices not allowed.")
  if (any(ordering != seq(1L,l))) slot(x, order_slot) <- ordering
  return(x)
}


#' Extract matrix
#'
#' @param x an ondisc_matrix.
#'
#' @return an in-memory version of x
extract_matrix <- function(x) {
  # First, determine which axis to index on; always index on the shorter axis.
  x_dim <- dim(x)
  if (x_dim[1] == 0 || x_dim[2] == 0) {
    stop("ondisc_matrix has 0 rows or columns; cannot extract sub-matrix into memory.")
  }
  index_on_cell <- x_dim[2] < x_dim[1]
  subset_vector <- get_subset_vector(x, index_on_cell)
  order_vector <- get_order_vector(x, index_on_cell)
  # if subset_vector is NA, initialize it as an integer vector (this will happen only rarely)
  if (identical(subset_vector, NA_integer_)) {
    subset_vector <- seq(1, if (index_on_cell) x_dim[2] else x_dim[1])
  }
  # if necessary, order the subset vector before passing to return_spMatrix_from_index.
  subset_vector_ordered <- identical(order_vector, NA_integer_)
  if (!subset_vector_ordered) {
    subset_vector <- subset_vector[order_vector]
  }
  # return submatrix, ordered (and subset) along main axis, and unsubset along secondary axis.
  out <- return_spMatrix_from_index(h5_file = x@h5_file, idx = subset_vector, col_idx = index_on_cell)
  # Two additional post-processing tasks:
  # First, reorder the rows/columns of the matrix (if necessary).
  if (!subset_vector_ordered) {
    inv_perm <- Matrix::invPerm(order_vector)
    out <- if (index_on_cell) out[,inv_perm,drop=FALSE] else out[inv_perm,,drop=FALSE]
  }
  # Second, if necessary, subset along the other axis
  second_subset <- get_subset_vector(x, !index_on_cell)
  if (!identical(second_subset, NA_integer_)) {
    out <- if (index_on_cell) out[second_subset,,drop=FALSE] else out[,second_subset,drop=FALSE]
  }
  return(out)
}


#' Convert an `ondisc_matrix` to a Seurat object
#'
#' Converts an ondisc_matrix to a Suerat object.
#'
#' @param x an ondisc_matrix
#' @param project (default "CreateSeuratObject") project name of the Seurat object
#' @param assay (default "RNA") name of the initial assay
#' @param cell_covariate_matrix (default none) cell-level covariate matrix
#'
#' @return a Seurat object
#' @export
#' @examples
#' # NOTE: You must create the HDF5 file "example.h5" to run this example.
#' # Navigate to the help file of "create_ondisc_matrix_from_10x_mtx"
#' # (via ?create_ondisc_matrix_from_10x_mtx), and execute the code in the example.
#' odm_fp <- system.file("extdata", "example.h5", package = "ondisc")
#' if (odm_fp != "") { # if required file exists, ...
#' odm <- ondisc_matrix(h5_file = odm_fp)
#' seurat_obj <- convert_to_seurat_object(odm)
#' }
convert_to_seurat_object <- function(x, project = "CreateSeuratObject", assay = "RNA", cell_covariate_matrix = NULL) {
  feature_ids <- get_feature_ids(x)
  cell_barcodes <- get_cell_barcodes(x)
  m <- extract_matrix(x)
  row.names(m) <- feature_ids
  colnames(m) <- cell_barcodes
  out <- Seurat::CreateSeuratObject(counts = m, project = project,
                                    assay = assay,
                                    meta.data = cell_covariate_matrix)
  return(out)
}
