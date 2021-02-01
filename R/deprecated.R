# depcrecated


# Summarize method, and more general apply/reduce methods
##########################################################

# setGeneric(name = "apply", def = function(X, MARGIN, FUN, ...) standardGeneric("apply"))


#' Apply an arbitrary function to all rows or columns
#' @param X an on_disc_matrix
#' @param MARGIN apply to rows (1) or columns (2)
#' @param FUN a function to apply
#' @param chunk_size number of rows (genes) or columns (cells) to load at a time; default 10000
#' @examples
#' # NOTE: You must create the HDF5 file "example.h5" to run this example.
#' # Navigate to the help file of "create_on_disc_matrix_from_10x_mtx"
#' # (via ?create_on_disc_matrix_from_10x_mtx), and execute the code in the example.
#' odm_fp <- system.file("extdata", "example.h5", package = "ondisc")
#' if (odm_fp != "") { # if required file exists, ...
#' odm <- on_disc_matrix(h5_file = odm_fp)
#' gene_means_and_sds <- apply(X = odm,
#' MARGIN = 1,
#' FUN = function(r) c(mean = mean(r), sd = sd(r)),
#' chunk_size = 500)
#' }
# setMethod(f = "apply",
#          signature = signature("on_disc_matrix"),
#          definition = function(X, MARGIN, FUN, chunk_size = 10000) {
#            closure <- function(chunk) apply(X = as.matrix(chunk), MARGIN = MARGIN, FUN = FUN)
#            on_disc_apply_across_chunks(x = X, col_apply = (MARGIN == 2), chunk_function = closure, chunk_size = chunk_size)
#          })


#' Combine results
#'
#' Takes a list of results (of the same type, i.e. vector or matrix) and combines them.
#'
#' @param res_list a list of results
#'
#' @return a vector or matrix containing the combined results
# combine_results <- function(res_list) {
#  if (is(res_list[[1]], "index")) { # if the output is a vector, c
#    out <- do.call(what = c, args = res_list)
#  }
#  if (is(res_list[[1]], "matrix")) { # if matrices, cbind
#    out <- do.call(what = cbind, args = res_list)
#  }
#  return(out)
# }


#' Print checkmark
#' prints a checkmark
#' @return NULL
# print_checkmark <- function() {
#  cat(crayon::green('  \u2713')); cat("\n")
# }


#' Get a chunks sequence
#'
#' Generates a sequence of cells/genes to load given a chunk size and the number of cells/genes.
#'
#' @param n_items number of items (genes or cells)
#' @param chunk_size size of chunk
#'
#' @return a list of indices v; for entry i, (v(i) + 1):(v(i+1)) provides a range.
# get_chunks_sequence <- function(n_items, chunk_size) {
#  chunk_size <- min(chunk_size, n_items)
#  init_seq <- seq(from = 0, to = n_items, by = chunk_size)
#  if (n_items %% chunk_size != 0) init_seq <- c(init_seq, n_items)
#  return(init_seq)
# }


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
# on_disc_apply_across_chunks <- function(x, col_apply, chunk_function, chunk_size) {
#  cat(paste0("Chunk size: ", chunk_size, "\n"))
#  sequence <- get_chunks_sequence(n_items = if (col_apply) ncol(x) else nrow(x), chunk_size = chunk_size)
#  res_list <- purrr::map(.x = 1:(length(sequence)-1), .f = function(i) {
#    curr_idx <- (sequence[i] + 1):(sequence[i + 1])
#    cat(paste0("Processing chunk ", crayon::blue(i), " of ", crayon::blue(length(sequence) - 1), ".\n"))
#    curr_m <- subset_by_gene_or_cell(x = x, idx = curr_idx, subset_on_cell = col_apply) %>% extract_matrix()
#    chunk_function(curr_m)
#  })
#  cat("All chunks processed."); print_checkmark()
#  # Finally, combine all the results
#  combine_results(res_list)
# }


#' Get essential gene- and cell-wise summary statistics
#'
#' Takes an on_disc_matrix that contains gene-by-cell expressions (x), and outputs the cell-specific and gene-specific covariate matrices.
#' @param x an on_disc_matrix
#' @param n_cells_to_process_at_time (optional; default 10000) number of cells to process at a time
#'
#' @return a list containing the cell-specific and gene-specific covariate matrices
#' @examples
#' # NOTE: You must create the HDF5 file "example.h5" to run this example.
#' # Navigate to the help file of "create_on_disc_matrix_from_10x_mtx"
#' # (via ?create_on_disc_matrix_from_10x_mtx), and execute the code in the example.
#' odm_fp <- system.file("extdata", "example.h5", package = "ondisc")
#' if (odm_fp != "") { # if required file exists, ...
#' odm <- on_disc_matrix(h5_file = odm_fp)
#' covariate_matrices <- summarize_expression_matrix(odm)
#' }
# summarize_expression_matrix <- function(x, n_cells_to_process_at_time = 10000) {
#  # Obtain the gene_names (constant across cell chunks)
#  gene_ids <- get_gene_names(x)
#  mt_genes <- grep(pattern = "^MT-", x = gene_ids)
#  sequence <- get_chunks_sequence(n_items = ncol(x), chunk_size = n_cells_to_process_at_time)
#  res_list <- purrr::map(.x = 1:(length(sequence)-1), .f = function(i) {
#    curr_idx <- (sequence[i] + 1):(sequence[i + 1])
#    cat(paste0("Processing chunk ", crayon::blue(i), " of ", crayon::blue(length(sequence) - 1), ".\n"))
#    curr_m <- subset_by_gene_or_cell(x = x, idx = curr_idx, subset_on_cell = TRUE) %>% extract_matrix()
#    # cell-wise total UMI count
#    cell_wise_umi_total <- Matrix::colSums(curr_m)
#    # cell-wise mitochondrial UMI count
#    if (length(mt_genes) >= 1) {
#      cell_wise_mito_umi_total <- Matrix::colSums(curr_m[mt_genes,,drop=FALSE])
#    } else {
#      cell_wise_mito_umi_total <- numeric(length = ncol(curr_m))
#    }
#    # gene-wise total UMI count
#    gene_wise_umi_total <-  Matrix::rowSums(curr_m)
#    # Convert curr_m to binary matrix
#    curr_m@x <- rep(1, length(curr_m@x))
#    # cell-wise n genes expressed
#    cell_wise_n_genes_exp <-  Matrix::colSums(curr_m)
#    # gene-wise n-cells expressed
#    gene_wise_n_cells_exp <-  Matrix::rowSums(curr_m)
#    # Return data in two matrices: a gene covaraite matrix and cell covariate matrix
#    gene_covariate_matrix <- cbind(gene_wise_umi_total, gene_wise_n_cells_exp)
#    cell_covariate_matrix <- cbind(cell_wise_umi_total, cell_wise_n_genes_exp, cell_wise_mito_umi_total)
#    list(gene_covariate_matrix = gene_covariate_matrix, cell_covariate_matrix = cell_covariate_matrix)
#  })
#  cat("All chunks processed."); print_checkmark()
#  # Finally, combine all the results: for genes, reduce over sum; for cells, do.call over rbind.
#  cell_covariate_matrix_temp <- do.call(rbind, purrr::map(.x = res_list, .f = function(i) i[["cell_covariate_matrix"]]))
#  p_mito <- (cell_covariate_matrix_temp[,"cell_wise_mito_umi_total"]/cell_covariate_matrix_temp[,"cell_wise_umi_total"]) * 100
#  cell_covariate_matrix <- cbind(cell_covariate_matrix_temp[,c("cell_wise_umi_total", "cell_wise_n_genes_exp")], p_mito = p_mito)
#  colnames(cell_covariate_matrix) <- c("total_umis", "n_genes_expressed", "percent_mito")
#  gene_cov_mat_dim <- dim(res_list[[1]][["gene_covariate_matrix"]])
#  gene_covariate_matrix <- Reduce(f = "+", x = purrr::map(.x = res_list, .f = function(i) i[["gene_covariate_matrix"]]), init = matrix(data = 0, nrow = gene_cov_mat_dim[1], ncol = gene_cov_mat_dim[2]))
#  colnames(gene_covariate_matrix) <- c("total_umis", "n_cells_expressed")
#  return(list(cell_covariate_matrix = cell_covariate_matrix, gene_covariate_matrix = gene_covariate_matrix))
# }
