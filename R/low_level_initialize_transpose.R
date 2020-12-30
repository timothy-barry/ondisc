#' Generate on disc_matrix_name
#'
#' Generates the name of an on_disc_matrix object given a directory. This function searches for files named on_disc_matrix_*.h5 in the specified directory. If none exists, it returns on_disc_matrix_1.h5. Else, it returns n_disc_matrix_*.h5 with a unique integer in place of *.
#'
#' @param on_disc_dir directory in which to store the on_disc_matrix.
#' @return a new name for an on_disc_matrix.
generate_on_disc_matrix_name <- function(on_disc_dir) {
  fs <- list.files(on_disc_dir)
  base_name <- "on_disc_matrix_"
  idxs <- grep(pattern = paste0(base_name, "[0-9]+.h5"), x = fs)
  if (length(idxs) == 0) {
    name <- paste0(base_name, "1.h5")
  } else {
    existing_names <- fs[idxs]
    ints_in_use <- gsub(pattern = paste0(base_name, "(\\d+).h5"), replacement = "\\1", x = existing_names) %>% as.integer()
    new_int <- max(ints_in_use) + 1
    name <- paste0(base_name, new_int, ".h5")
  }
  return(paste0(on_disc_dir, "/", name))
}

#' Create h5 file on disk
#'
#' @param on_disc_dir directory in which to store the on_disc object
#' @param n_genes number of genes in dataset
#' @param n_cells number of cells in dataset
#' @param n_data_points total number of datapoints
#' @param cell_barcodes cell barcodes
#' @param gene_ids the gene ids
#' @param gene_names the gene names
#' @param file_name name if the file to create; if null, generate_on_disc_matrix_name is called.
#' @return location of h5 file on disk.
create_h5_file_on_disk <- function(on_disc_dir, n_genes, n_cells, n_data_points, cell_barcodes, gene_ids, gene_names, file_name) {
  cat("Initializing the on_disc_matrix...")
  on_disc_location <- if (is.null(file_name)) generate_on_disc_matrix_name(on_disc_dir) else paste0(on_disc_dir, "/", file_name, ".h5")
  # Remove the -1 tags from the cell barcodes
  cell_barcodes <- gsub(pattern = '*-1', replacement = "", x = cell_barcodes)
  # Initialize the group structure.
  rhdf5::h5createFile(on_disc_location) %>% invisible()
  rhdf5::h5createGroup(file = on_disc_location, group = "metadata") %>% invisible()
  rhdf5::h5createGroup(file = on_disc_location, group = "compressed_sparse_column") %>% invisible()
  rhdf5::h5createGroup(file = on_disc_location, group = "compressed_sparse_row") %>% invisible()
  # Write the metadata vectors.
  rhdf5::h5write(cell_barcodes, on_disc_location, "metadata/cell_barcodes")
  rhdf5::h5write(gene_ids, on_disc_location, "metadata/gene_ids")
  rhdf5::h5write(gene_names, on_disc_location, "metadata/gene_names")
  dimension <- c(n_genes, n_cells)
  rhdf5::h5write(dimension, on_disc_location, "metadata/dimension")
  # Initialize the data vectors in the compressed sparse column and row groups.
  rhdf5::h5createDataset(file = on_disc_location, dataset = "compressed_sparse_column/data", dims = c(n_data_points, 2), storage.mode = "integer", level = 0, chunk = c(100, 1)) %>% invisible()
  rhdf5::h5createDataset(file = on_disc_location, dataset = "compressed_sparse_column/p", dims = n_cells + 1, storage.mode = "integer", level = 0, chunk = 10) %>% invisible()

  rhdf5::h5createDataset(file = on_disc_location, dataset = "compressed_sparse_row/data", dims = c(n_data_points, 2), storage.mode = "integer", level = 0, chunk = c(100, 1)) %>% invisible()
  rhdf5::h5createDataset(file = on_disc_location, dataset = "compressed_sparse_row/p", dims = n_genes + 1, storage.mode = "integer", level = 0, chunk = 10) %>% invisible()
  print_checkmark()
  # Return name of h5 file
  return(on_disc_location)
}


#' Increment indices
#'
#' A simple helper function to increment the indices of a vector. (Note: re-write in C++)
#'
#' @param old_indices A vector
#' @param new_indices the indices (optionally repeated) to increment by 1 in the vector
#'
#' @return the old_indices vector with appropriately updated values.
increment_indices <- function(old_indices, new_indices) {
  for (curr_idx in new_indices) {
    old_indices[curr_idx] <- old_indices[curr_idx] + 1
  }
  return(old_indices)
}


#' Convert mtx file to compressed sparse column
#'
#' Indexing note: the input (.mtx) uses 1-based indexing for cell and gene. This function saves the (unchanged) expression data; the gene indexes, which are converted to 0-based indexing; and the row and column pointers, which are saved using 1-based indexing.
#'
#' @param mtx_fp file path to mtx file
#' @param h5_loc file path to (already-initialized) h5 file
#' @param n_rows_with_comments number of rows with comments in the mtx file; the first n_rows_with_comments + 1 rows will be skipped.
#' @param chunk_size chunk size for the read_delim_chunked function; default 1e5
#'
#' @return NULL.
convert_mtx_to_csc <- function(mtx_fp, h5_loc, n_rows_with_comments, chunk_size = 1e5) {
  # increment n_rows_with_comments to avoid the row containing the number of genes, cells, and data points.
  n_rows_with_comments <- n_rows_with_comments + 1
  # recall the number of genes and cells
  dimension <- rhdf5::h5read(h5_loc, "metadata/dimension")
  n_genes <- dimension[1]
  n_cells <- dimension[2]
  # Initialize the n_data_points_gene and n_data_points_cell accumulator vectors.
  acc_init <- list(n_data_points_gene = integer(length = n_genes), n_data_points_cell = integer(length = n_cells))
  # Define the closure to be called by read_delim_chunked
  f <- function(x, pos, acc) {
    # First, append the data to the h5 file.
    n_elements <- nrow(x)
    curr_index <- pos:(pos + n_elements - 1)
    rhdf5::h5write(obj = dplyr::select(x, gene_idx, expression) %>% dplyr::mutate(gene_idx = gene_idx - 1) %>% as.matrix(), file = h5_loc, name = "compressed_sparse_column/data", index = list(curr_index, NULL))
   # Next, increment the gene index and cell index counters.
    n_data_points_gene <- increment_indices(old_indices = acc$n_data_points_gene, new_indices = x$gene_idx)
    n_data_points_cell <- increment_indices(old_indices = acc$n_data_points_cell, new_indices = x$cell_idx)
    return(list(n_data_points_gene = n_data_points_gene, n_data_points_cell = n_data_points_cell))
  }
  data_indices <- readr::read_delim_chunked(file = mtx_fp, chunk_size = chunk_size, skip = n_rows_with_comments, callback = readr::AccumulateCallback$new(f, acc = acc_init), delim = " ", col_names = c("gene_idx", "cell_idx", "expression"), progress = TRUE, col_types = "iii")
  col_ptr <- c(1, cumsum(data_indices$n_data_points_cell) + 1)
  row_ptr <- c(1, cumsum(data_indices$n_data_points_gene) + 1)
  rhdf5::h5write(obj = col_ptr, file = h5_loc, name = "compressed_sparse_column/p")
  rhdf5::h5write(obj = row_ptr, file = h5_loc, name = "compressed_sparse_row/p")
  return(invisible())
}


#' Transpose on-disc CSC matrix.
#'
#' Converts a CSC matrix to a CSR matrix. The CSC matrix is assumed to be stored in an h5 file with fields compressed_sparse_column/p, compressed_sparse_column/data, and compressed_sparse_row/p. The function writes compressed_sparse_row/data. Indexing note: The row and column pointers are assumed to be stored using 1-based indexing. The gene indexes, as stored in compressed_sparse_column/data, use 0-based indexing. The cell indexes to be saved to compressed_sparse_row/data will be saved using 0-based indexing.
#'
#' @param h5_file an h5 file on disk.
#' @param cell_chunk_size Load chunks of the data containing this many cells.
#'
#' @return invisible
transpose_on_disc_csc_matrix <- function(h5_file, cell_chunk_size = 5000) {
  # Load pointers and basic metadata
  row_ptr <- rhdf5::h5read(file = h5_file, name = "/compressed_sparse_row/p")
  col_ptr <- rhdf5::h5read(file = h5_file, name = "/compressed_sparse_column/p")
  dim <- rhdf5::h5read(file = h5_file, name = "/metadata/dimension")
  n_cells <- dim[2]; n_genes <- dim[1]

  # create the cell sequence
  cell_chunk_sequence <- get_chunks_sequence(n_cells, cell_chunk_size)
  n_chunks <- length(cell_chunk_sequence) - 1
  single_run <- n_chunks == 1

  # Determine number of entries per row and column
  n_entries_per_col <- col_ptr %>% diff()

  # If single_run is true, simply load all data into memory, sort by gene id cell id, and save the result.
  if (single_run) {
    cat("Loading data from disk.\n")
    dt <- data.table::data.table(rep(x = 0:(n_cells - 1), times = n_entries_per_col),
                                 rhdf5::h5read(file = h5_file, name = "/compressed_sparse_column/data"))
    colnames(dt) <- c("cell_idx", "gene_idx", "expression")
    data.table::setorderv(dt, c("gene_idx", "cell_idx"))
    # save cell index-expression matrix
    rhdf5::h5write(file = h5_file, obj = dt[, c("cell_idx", "expression")] %>% as.matrix(), name = "/compressed_sparse_row/data", index = list(NULL, NULL))
    return(invisible())
  }

  # Perform the oom matrix transpose algorithm: Load the chunked data, sort by gene id, increment the row pointer, and store the results.
  for (chunk in 1:n_chunks) {
    cell_start <- cell_chunk_sequence[chunk] + 1
    cell_end <- cell_chunk_sequence[chunk + 1]
    # determine which part of CSC matrix to load
    col_ptr_range <- col_ptr[c(cell_start, cell_end + 1)] + c(0, -1)
    # Create triplet matrix
    cat("Loading chunk of data from disk.\n")
    dt <- data.table::data.table(rep(x = (cell_start - 1):(cell_end-1), times = n_entries_per_col[cell_start:cell_end]),
                                 rhdf5::h5read(file = h5_file, name = "/compressed_sparse_column/data", index = list(col_ptr_range[1]:col_ptr_range[2], NULL)))
    cat(paste("Data chunk", chunk, "of", n_chunks, "loaded.\n"))
    colnames(dt) <- c("cell_idx", "gene_idx", "expression")
    # arrange the data table by gene_idx
    data.table::setorderv(dt, c("gene_idx", "cell_idx"))
    # determine number of expressed genes per gene id
    n_entries_per_gene <- tabulate(dt$gene_idx + 1)
    nonzero_genes <- n_entries_per_gene >= 1
    n_entries_per_gene_nonzero <- n_entries_per_gene[nonzero_genes]
    gene_id_nonzero <- (1:n_genes)[nonzero_genes]
    # Determine the row indices of the data to store in CSR format
    start_end_locs <- purrr::map(.x = 1:length(n_entries_per_gene_nonzero), .f = function(i) {
      start_loc <- row_ptr[[gene_id_nonzero[i]]]
      end_loc <- start_loc + n_entries_per_gene_nonzero[i] - 1
      c(start_loc = start_loc, end_loc = end_loc)
    })
    # Obtain the cumulative sum of the gene counts
    gene_cum_sum <- c(0, cumsum(n_entries_per_gene_nonzero))
    # Iterate over the genes, laying down the data in the appropriate position in the CSR matrix.
    for (i in 1:length(n_entries_per_gene_nonzero)) {
      if (i %% 1000 == 0) cat(paste0("Writing gene ", i, " of ", length(start_end_locs), " to disk.\n"))
      inds <- (gene_cum_sum[i] + 1):gene_cum_sum[i+1]
      # Lay down the cell indexes and expression
      rhdf5::h5write(file = h5_file, obj = dt[inds, c("cell_idx", "expression")] %>% as.matrix(), name = "/compressed_sparse_row/data", index = list(start_end_locs[[i]][1]:start_end_locs[[i]][2], NULL))
      # increment the row pointer
      row_ptr[gene_id_nonzero[i]] <- row_ptr[gene_id_nonzero[i]] + n_entries_per_gene_nonzero[i]
    }
  }
  return(invisible())
}
