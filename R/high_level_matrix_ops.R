#' Get on_disc subset vector
#'
#' Returns the subset vector (either gene or cell) for an on_disc_matrix object.
#'
#' @param x the on_disc_matrix.
#' @param cell_idx (boolean) TRUE = cell, FALSE = gene.
#'
#' @return the requested subset indexes
get_subset_vector <- function(x, cell_idx) {
  idx_slot <- paste0(if (cell_idx) "cell" else "gene", "_subset")
  return(slot(x, idx_slot))
}


#' Get dim
#'
#' Returns the dimension of an on_disc matrix.
#'
#' @param x an on_disc_matrix object
#'
#' @return the dimension
get_dim <- function(x) {
  cell_unsubset <- identical(x@cell_subset, NA_integer_)
  gene_unsubset <- identical(x@gene_subset, NA_integer_)

  if (cell_unsubset || gene_unsubset) dim_h5 <- rhdf5::h5read(file = x@h5_file, name = "/metadata/dimension") %>% as.integer()
  n_row_out <- if (gene_unsubset) dim_h5[1] else length(x@gene_subset)
  n_col_out <- if (cell_unsubset) dim_h5[2] else length(x@cell_subset)
  return(c(n_row_out, n_col_out))
}


#' Get names helper function
#'
#' @param x an on_disc_matrix
#' @param name_to_get the names to extract from x (one of cell_barcodes, gene_ids, gene_names)
#'
#' @return a character vector containing the requested names
get_names <- function(x, name_to_get) {
  names_out <- rhdf5::h5read(file = x@h5_file, name = paste0("metadata/", name_to_get)) %>% as.character()
  idx <- get_subset_vector(x, name_to_get == "cell_barcodes")
  if (identical(idx, NA_integer_)) names_out else names_out[idx]
}


#' Get names
#'
#' Functions to extract cell barcodes, gene ids, and gene names from an on_disc_matrix.
#'
#' @param x An on_disc_matrix.
#' @return The cell barcodes of this on-disc matrix.
#' @export
get_cell_barcodes <- function(x) {
  get_names(x, "cell_barcodes")
}


#' @rdname get_cell_barcodes
#' @param x An on_disc_matrix.
#' @export
#' @return The cell gene ids of this on-disc matrix.
get_gene_ids <- function(x) {
  get_names(x, "gene_ids")
}


#' @rdname get_cell_barcodes
#' @param x An on_disc_matrix.
#' @export
#' @return The gene names of this on-disc matrix.
get_gene_names <- function(x) {
  get_names(x, "gene_names")
}


#' Subset by gene or cell
#'
#' Subsets an on_disc_matrix by either gene or cell, as specified.
#'
#' @param x an on_disc_matrix
#' @param idx a set of integer, logical, or character indexes.
#' @param subset_on_cell boolean indicating whether to subset on cell (TRUE) or gene (FALSE)
#'
#' @return a subset on_disc_matrix.
subset_by_gene_or_cell <- function(x, idx, subset_on_cell) {
  n_elements <- if (subset_on_cell) ncol(x) else nrow(x)
  subset_slot <- paste0(if (subset_on_cell) "cell" else "gene", "_subset")
  if (identical(slot(x, subset_slot), NA_integer_)) {
    slot(x, subset_slot) <- 1:n_elements
  }
  # Ensure input not too long
  if (length(idx) > n_elements) stop("Index vector too long.")
  # perform subset based on class of idx
  if (class(idx) %in% c("logical", "numeric", "integer")) {
    # Check if index is out of bounds
    if (max(idx) > length(slot(x, subset_slot))) stop("Index out of bounds.")
    slot(x, subset_slot) <- slot(x, subset_slot)[idx]
  }
  if (class(idx) == "character") {
    all_ids <- rhdf5::h5read(file = x@h5_file, name = paste0("metadata/", if (subset_on_cell) "cell_barcodes" else "gene_ids"))
    curr_ids <- all_ids[slot(x, subset_slot)]
    within_cur_ids_idxs <- match(x = idx, table = curr_ids)
    if (any(is.na(within_cur_ids_idxs))) stop(paste(if(subset_on_cell) "Cell barcode" else "Gene id","not present in data."))
    slot(x, subset_slot) <- slot(x, subset_slot)[within_cur_ids_idxs]
  }
  # Finally, verify that the indexes are not duplicated
  tmp <- rle(sort(slot(x, subset_slot), method = "radix"))
  if (any(tmp$lengths >= 2)) stop("Duplicate indices not allowed.")
  return(x)
}


#' Extract matrix
#'
#' @param x an on_disc matrix
#'
#' @return a sparse Matrix object corresponding to the in-memory version of x.
extract_matrix <- function(x) {
  # First, determine which axis to index on; always index on the shorter axis.
  x_dim <- dim(x)
  if (x_dim[1] == 0 || x_dim[2] == 0) stop("On_disc_matrix has 0 rows or columns; cannot extract sub-matrix into memory.")
  index_on_cell <- x_dim[2] < x_dim[1]
  # In rare cases, the subset vector may be NA; for now, just initialze the subset vector of the (deep-copied) on_disc_matrix to 1:nrow or 1:ncol. Later, maybe optimize by passing NULL to rhdf5 to load entire dataset.
  curr_subset <- paste0(if (index_on_cell) "cell" else "gene", "_subset")
  if (identical(slot(x, curr_subset), NA_integer_)) slot(x, curr_subset) <- if (index_on_cell) 1:x_dim[2] else 1:x_dim[1]
  # Next, pull data from memory according to the correct matrix axis.
  out <- return_spMatrix_from_index(h5_file = x@h5_file, idx = if (index_on_cell) x@cell_subset else x@gene_subset, col_idx = index_on_cell)
  # Next, subset the sparse matrix according to the other dimension.
  second_subset <- get_subset_vector(x, !index_on_cell)
  if (!identical(second_subset, NA_integer_)) { # If we have to perform a subset,
    if (index_on_cell) { # then subset by gene by we extracted by cell.
      out <- out[second_subset,,drop=FALSE]
    } else { # Or subset by cell if we extracted by gene.
      out <- out[,second_subset,drop=FALSE]
    }
  }
  return(out)
}


#' On disc apply across chunks
#'
#' Applies a function across chunks of an on_disc_matrix.
#'
#' @param x an on_disc matrix
#' @param col_apply (boolean) apply across columns (TRUE) or rows (FALSE)
#' @param chunk_function function to apply to each chunk
#' @param chunk_size size of the chunks to load
#'
#' @return a vector or matrix containing the results of applying the function to the on_disc matrix.
on_disc_apply_across_chunks <- function(x, col_apply, chunk_function, chunk_size) {
  cat(paste0("Chunk size: ", chunk_size, "\n"))
  sequence <- get_chunks_sequence(n_items = if (col_apply) ncol(x) else nrow(x), chunk_size = chunk_size)
  res_list <- purrr::map(.x = 1:(length(sequence)-1), .f = function(i) {
    curr_idx <- (sequence[i] + 1):(sequence[i + 1])
    cat(paste0("Processing chunk ", crayon::blue(i), " of ", crayon::blue(length(sequence) - 1), ".\n"))
    curr_m <- subset_by_gene_or_cell(x = x, idx = curr_idx, subset_on_cell = col_apply) %>% extract_matrix()
    chunk_function(curr_m)
  })
  cat("All chunks processed."); print_checkmark()
  # Finally, combine all the results
  combine_results(res_list)
}


#' Summarize expression matrix
#'
#' Takes an on_disc_matrix that contains gene-by-cell expressions (x), and outputs the cell-specific and gene-specific covariate matrices.
#'
#' @param x an on_disc_matrix
#' @param chunk_size number of cells to process at a time
#'
#' @return a list containing the cell-specific and gene-specific covariate matrices
#' @export
#' @examples
#' exp_mat_loc <- system.file("extdata", "on_disc_matrix_1.h5", package = "ondisc")
#' if (exp_mat_loc != "") {
#' x <- on_disc_matrix(h5_file = exp_mat_loc)
#' covariate_matrices <- summarize_expression_matrix(x)
#' }
summarize_expression_matrix <- function(x, chunk_size = 4000) {
  # Obtain the gene_names (constant across cell chunks)
  gene_ids <- get_gene_names(x)
  mt_genes <- grep(pattern = "^MT-", x = gene_ids)
  sequence <- get_chunks_sequence(n_items = ncol(x), chunk_size = chunk_size)
  res_list <- purrr::map(.x = 1:(length(sequence)-1), .f = function(i) {
    curr_idx <- (sequence[i] + 1):(sequence[i + 1])
    cat(paste0("Processing chunk ", crayon::blue(i), " of ", crayon::blue(length(sequence) - 1), ".\n"))
    curr_m <- subset_by_gene_or_cell(x = x, idx = curr_idx, subset_on_cell = TRUE) %>% extract_matrix()
    # cell-wise total UMI count
    cell_wise_umi_total <- Matrix::colSums(curr_m)
    # cell-wise mitochondrial UMI count
    if (length(mt_genes) >= 1) {
      cell_wise_mito_umi_total <-  Matrix::colSums(curr_m[mt_genes,])
    } else {
      cell_wise_mito_umi_total <- numeric(length = ncol(x))
    }
    # gene-wise total UMI count
    gene_wise_umi_total <-  Matrix::rowSums(curr_m)
    # Convert curr_m to binary matrix
    curr_m@x <- rep(1, length(curr_m@x))
    # cell-wise n genes expressed
    cell_wise_n_genes_exp <-  Matrix::colSums(curr_m)
    # gene-wise n-cells expressed
    gene_wise_n_cells_exp <-  Matrix::rowSums(curr_m)
    # Return data in two matrices: a gene covaraite matrix and cell covariate matrix
    gene_covariate_matrix <- cbind(gene_wise_umi_total, gene_wise_n_cells_exp)
    cell_covariate_matrix <- cbind(cell_wise_umi_total, cell_wise_n_genes_exp, cell_wise_mito_umi_total)
    list(gene_covariate_matrix = gene_covariate_matrix, cell_covariate_matrix = cell_covariate_matrix)
  })
  cat("All chunks processed."); print_checkmark()
  # Finally, combine all the results: for genes, reduce over sum; for cells, do.call over rbind.
  cell_covariate_matrix_temp <- do.call(rbind, purrr::map(.x = res_list, .f = function(i) i[["cell_covariate_matrix"]]))
  p_mito <- (cell_covariate_matrix_temp[,"cell_wise_mito_umi_total"]/cell_covariate_matrix_temp[,"cell_wise_umi_total"]) * 100
  cell_covariate_matrix <- cbind(cell_covariate_matrix_temp[,c("cell_wise_umi_total", "cell_wise_n_genes_exp")], p_mito = p_mito)
  colnames(cell_covariate_matrix) <- c("total_umis", "n_genes_expressed", "percent_mito")
  gene_cov_mat_dim <- dim(res_list[[1]][["gene_covariate_matrix"]])
  gene_covariate_matrix <- Reduce(f = "+", x = purrr::map(.x = res_list, .f = function(i) i[["gene_covariate_matrix"]]), init = matrix(data = 0, nrow = gene_cov_mat_dim[1], ncol = gene_cov_mat_dim[2]))
  colnames(gene_covariate_matrix) <- c("total_umis", "n_cells_expressed")
  return(list(cell_covariate_matrix = cell_covariate_matrix, gene_covariate_matrix = gene_covariate_matrix))
}
