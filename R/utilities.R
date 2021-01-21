#' Print checkmark
#' prints a checkmark
#' @return NULL
print_checkmark <- function() {
  cat(crayon::green('  \u2713')); cat("\n")
}


#' Get a chunks sequence
#'
#' Generates a sequence of cells/genes to load given a chunk size and the number of cells/genes.
#'
#' @param n_items number of items (genes or cells)
#' @param chunk_size size of chunk
#'
#' @return a list of indices v; for entry i, (v(i) + 1):(v(i+1)) provides a range.
get_chunks_sequence <- function(n_items, chunk_size) {
  chunk_size <- min(chunk_size, n_items)
  init_seq <- seq(from = 0, to = n_items, by = chunk_size)
  if (n_items %% chunk_size != 0) init_seq <- c(init_seq, n_items)
  return(init_seq)
}


#' Combine results
#'
#' Takes a list of results (of the same type, i.e. vector or matrix) and combines them.
#'
#' @param res_list a list of results
#'
#' @return a vector or matrix containing the combined results
combine_results <- function(res_list) {
  if (is(res_list[[1]], "index")) { # if the output is a vector, c
    out <- do.call(what = c, args = res_list)
  }
  if (is(res_list[[1]], "matrix")) { # if matrices, cbind
    out <- do.call(what = cbind, args = res_list)
  }
  return(out)
}


#' Get n rows with comments
#'
#' Returns the number of rows with comments in an mtx file.
#'
#' @param mtx_fp a file path to an mtx file
#'
#' @return the number of rows with comments at the top of the file
get_n_rows_with_comments_mtx <- function(mtx_fp) {
  n_rows_with_comments <- 0
  repeat {
    curr_row <- utils::read.table(mtx_fp, nrows = 1, skip = n_rows_with_comments, header = FALSE, sep = "\n") %>% dplyr::pull()
    is_comment <- substr(curr_row, start = 1, stop = 1) == "%"
    if (!is_comment) {
      break()
    } else {
      n_rows_with_comments <- n_rows_with_comments + 1
    }
  }
  return(n_rows_with_comments)
}


#' Get mtx metadata
#'
#' @param mtx_fp filepath to the mtx file
#' @param n_rows_with_comments number of rows with comments (at top of file)
#'
#' @return a list containing (i) n_genes, (ii) n_cells, (iii) the sparsity (i.e., fraction of entries that are zero), (iv) (TRUE/FALSE) matrix is logical
get_mtx_metadata <- function(mtx_fp, n_rows_with_comments) {
  metadata <- utils::read.table(file = mtx_fp, nrows = 1, skip = n_rows_with_comments, header = FALSE, sep = " ", colClasses = c("integer", "integer", "integer"))
  n_genes <- metadata %>% dplyr::pull(1)
  n_cells <- metadata %>% dplyr::pull(2)
  n_data_points <- metadata %>% dplyr::pull(3)
  sparsity <- 1 - n_entries / (n_genes * n_cells)
  first_row <- utils::read.table(file = mtx_fp, nrows = 1, skip = n_rows_with_comments + 1, header = FALSE, sep = " ", colClasses = c("integer", "integer", "integer"))
  is_logical <- ncol(first_row) == 2
  return(list(n_genes = n_genes, n_cells = n_cells, n_data_points = n_data_points, sparsity = sparsity, is_logical = is_logical))
}


#' n GB to entries
#'
#' @param n_gb number of gigabytes to process per chunk
#' @param logical_mtx number of
#'
#' @return
#' @export
#'
#' @examples
n_gb_to_n_entries <- function(n_gb, logical_mtx) {
  multiplicative_factor <- 1e9 * (if (logical_mtx) 8.0 else 12.0)
  return (multiplicative_factor * n_gb)
}
