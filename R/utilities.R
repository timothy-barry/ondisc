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
