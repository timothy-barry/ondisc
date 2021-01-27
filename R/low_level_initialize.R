#' Initialize h5 file on-disk
#'
#' Initialize the on-disk portion on an ondisc_matrix.
#'
#' @param h5_fp file path to the .h5 file to be initialized
#' @param mtx_metadata metadata of the .mtx file
#' @param features_metadata metadata of the features.tsv file
#' @param barcodes_fp file path to the barcodes.tsv file
#' @param features_fp file path to the features.tsv file
#'
#' @return NULL
initialize_h5_file_on_disk <- function(h5_fp, mtx_metadata, features_metadata, barcodes_fp, features_fp) {
  # Create the .h5 file
  rhdf5::h5createFile(h5_fp) %>% invisible()
  # Write metadata
  cell_barcodes <- dplyr::pull(readr::read_tsv(file = barcodes_fp, col_names = FALSE, col_types = "c"))
  rhdf5::h5write(cell_barcodes, h5_fp, "cell_barcodes")
  feature_ids <- read_given_column_of_tsv(col_idx = 1, n_cols = features_metadata$n_cols, tsv_file = features_fp)
  rhdf5::h5write(feature_ids, h5_fp, "feature_ids")
  if (features_metadata$feature_names) {
    feature_names <- read_given_column_of_tsv(col_idx = 2, n_cols = features_metadata$n_cols, tsv_file = features_fp)
    rhdf5::h5write(feature_names, h5_fp, "feature_names")
  }
  rhdf5::h5write(c(mtx_metadata$n_features, mtx_metadata$n_cells), h5_fp, "dimension")
  # Initialize CSC
  rhdf5::h5createDataset(file = h5_fp, dataset = "col_ptr", dims = mtx_metadata$n_cells + 1, storage.mode = "integer", level = 0, chunk = 10) %>% invisible()
  rhdf5::h5createDataset(file = h5_fp, dataset = "row_idxs", dims = mtx_metadata$n_data_points, storage.mode = "integer", level = 0, chunk = 100) %>% invisible()
  if (!mtx_metadata$is_logical) {
  rhdf5::h5createDataset(file = h5_fp, dataset = "data_csc", dims = mtx_metadata$n_data_points, storage.mode = "integer", level = 0, chunk = 100) %>% invisible()
  }
  # Initialize CSR
  rhdf5::h5createDataset(file = h5_fp, dataset = "row_ptr", dims = mtx_metadata$n_features + 1, storage.mode = "integer", level = 0, chunk = 10) %>% invisible()
  rhdf5::h5createDataset(file = h5_fp, dataset = "col_idxs", dims = mtx_metadata$n_data_points, storage.mode = "integer", level = 0, chunk = 100) %>% invisible()
  if (!mtx_metadata$is_logical) {
    rhdf5::h5createDataset(file = h5_fp, dataset = "data_csr", dims = mtx_metadata$n_data_points, storage.mode = "integer", level = 0, chunk = 100) %>% invisible()
  }
  return(invisible())
}


#' Initialize accumulator
#'
#' @param terminal_symbol the terminal symbol
#' @param bag_of_variables the bag of variables
#'
#' @return An initialized vector of the correct type and length
initialize_accumulator <- function(terminal_symbol, bag_of_variables) {
  info <- get_terminal_acc_and_args(terminal_symbol)
  info$acc_constructor(bag_of_variables[[info$acc_length]])
}


#' Get accumulator funct and arg list
#'
#' @param terminal_symbol a terminal symbol
#'
#' @return a list containing (i) the name of and (ii) the arguments to pass to the accumulator function.
get_accumulator_funct_arg_list <- function(terminal_symbol) {
  info <- get_terminal_acc_and_args(terminal_symbol)
  list(acc_funct = info$acc_funct, acc_args = info$acc_args)
}


#' Run subtask 2a
#'
#' Runs subtask a of part 2 of mtx algorithm: updates the accumulator of each terminal symbol.
#'
#' @param x a data.table (passed by ref)
#' @param bag_of_variables the bag_of_variables (also passed by ref)
#' @param acc the accumulator list
#'
#' @return NULL
run_subtask_2a <- function(x, bag_of_variables, acc) {
  arguments <- arguments_enum()
  bag_of_variables[[arguments$cell_idxs]] <- x$cell_idxs
  bag_of_variables[[arguments$feature_idxs]] <- x$feature_idxs
  bag_of_variables[[arguments$umi_counts]] <- x$umi_counts
  for (i in 1:length(acc)) {
    # Add current acc_vect to bag of variables
    bag_of_variables[[arguments$acc_vect]] <- acc[[i]]
    # Obtain relevant arguments via mget
    curr_args <- mget(x = c(arguments$acc_vect, terminal_functs_args[[i]]$acc_args),
                      envir = bag_of_variables) %>% setNames(NULL)
    # Call relevant function
    do.call(what = terminal_functs_args[[i]]$acc_funct, args = curr_args)
  }
  return(invisible())
}


#' Run mtx algo step 2
#'
#' This function runs step 2 of the core mtx algorithm. It (a) computes the terminal symbols, (b) writes to the CSC matrix, and (b) sorts the data by feature_idx, then writes to the CSR matrix.
#'
#' @param h5_fp full path to the h5 file on-disk
#' @param n_elem_per_chunk number of elements to load per chunk
#' @param bag_of_variables the bag of variables containing the variables to pass to the accumulator functions
#' @param initial_accumulators list of starting accumulators
#' @param terminal_functs_args list of accumulator function names and arguments
#'
#' @return
run_mtx_algo_step_2 <- function(h5_fp, mtx_fp, is_logical, bag_of_variables, n_elem_per_chunk, n_rows_to_skip, initial_accumulators, terminal_functs_args) {
  # Define closure to be called by readr::read_delim_chunked
  closure <- function(x, pos, acc) {
    # example chunk: x <- read.table(file = mtx_fp, header = FALSE, sep = " ", col.names = c("feature_idxs", "cell_idxs", if (is_logical) NULL else "umi_counts"), skip = n_rows_to_skip, colClasses = rep("integer", if (is_logical) 2 else 3), nrows = 1000); pos <- 1000
    data.table::setDT(x)
    decrement_idxs(x$feature_idxs)
    decrement_idxs(x$cell_idxs)

    # subtask a: loop through accumulators and update
    run_subtask_2a(x, bag_of_variables, acc)

    # subtask b:
    # subtask c:
    return(acc)
  }

  terminals <- readr::read_delim_chunked(file = mtx_fp,
                            chunk_size = n_elem_per_chunk,
                            skip = n_rows_to_skip,
                            callback = readr::AccumulateCallback$new(closure, acc = initial_accumulators),
                            delim = " ",
                            col_names = c(arguments$feature_idxs, arguments$cell_idxs, if (is_logical) NULL else arguments$umi_counts ),
                            progress = TRUE,
                            col_types = if (is_logical) "ii" else "iii")
  return(terminals)
}


run_core_mtx_algo <- function(h5_fp, mtx_fp, is_logical, covariates, bag_of_variables, n_elem_per_chunk, n_rows_to_skip) {
  grammar <- initialize_grammar()

  # Determine which terminal symbols to compute
  terminal_symbols <- lapply(unlist(covariates),
                             get_terminals_for_covariate, grammar = grammar) %>% unlist() %>% unique()

  # Obtain the starting accumulator, as well as the accumulator function name and args, for each terminal
  initial_accumulators <- lapply(terminal_symbols, initialize_accumulator, bag_of_variables = bag_of_variables)
  terminal_functs_args <- lapply(terminal_symbols, get_accumulator_funct_arg_list)

  # Run step 1 of core algorithm

  # Run step 2 of core algorithm
  terminals <- run_mtx_algo_step_2(h5_fp, mtx_fp, is_logical, bag_of_variables, n_elem_per_chunk,
                                   n_rows_to_skip, initial_accumulators, terminal_functs_args)

  # compute the covariate matrices

}
