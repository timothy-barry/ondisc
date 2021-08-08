#' Initialize h5 file on-disk
#'
#' Initialize the on-disk portion on an ondisc_matrix.
#'
#' @param h5_fp file path to the .h5 file to be initialized
#' @param mtx_metadata metadata of the .mtx file
#' @param features_metadata metadata of the features.tsv file
#' @param barcodes file path to the barcodes.tsv file OR a barcodes R object
#' @param features file path to the features.tsv file OR a features R object
#' @param file_path bool argument to determine whether the input are file path or R object. TRUE for file path, FALSE for R object.
#'
#' @return NULL
#' @noRd
initialize_h5_file_on_disk <- function(h5_fp, mtx_metadata, features_metadata, barcodes, features, progress, file_path) {
  # Create the .h5 file
  status <- rhdf5::h5createFile(h5_fp)
  if(!status)
    stop(sprintf("Creating %s failed", h5_fp))

  # Write metadata
  cell_barcodes <- vector(mode = "character")
  if (file_path == TRUE) {
    if (length(barcodes) == 1) {
      cell_barcodes <- dplyr::pull(readr::read_tsv(file = barcodes, col_names = FALSE, col_types = "c", progress = progress))
    } else {
      # handle list of barcodes input
      for (i in 1:length(barcodes)) {
        single_cell_barcodes <- dplyr::pull(readr::read_tsv(file = barcodes[i], col_names = FALSE, col_types = "c", progress = progress))
        # concatenate cell_barcodes in each file to a single barcodes vector
        cell_barcodes <- c(cell_barcodes, single_cell_barcodes)
      }
    }
  } else {
    cell_barcodes <- barcodes
  }
  rhdf5::h5write(cell_barcodes, h5_fp, "cell_barcodes")

  if (file_path == TRUE) {
    feature_ids <- read_given_column_of_tsv(col_idx = 1, n_cols = features_metadata$n_cols, tsv_file = features, progress = progress)
  } else {
    feature_ids <- dplyr::pull(features, 1)
  }
  rhdf5::h5write(feature_ids, h5_fp, "feature_ids")

  if (features_metadata$feature_names) {
    if (file_path == TRUE) {
      feature_names <- read_given_column_of_tsv(col_idx = 2, n_cols = features_metadata$n_cols, tsv_file = features, progress = progress)
    } else {
      feature_names <- dplyr::pull(features, 2)
    }
    rhdf5::h5write(feature_names, h5_fp, "feature_names")
  }
  rhdf5::h5write(c(mtx_metadata$n_features, mtx_metadata$n_cells), h5_fp, "dimension")
  rhdf5::h5write(mtx_metadata$is_logical, h5_fp, "logical_mat")

  # Initialize CSC
  rhdf5::h5createDataset(file = h5_fp, dataset = "cell_ptr", dims = mtx_metadata$n_cells + 1, storage.mode = "integer", level = 0L, chunk = min(mtx_metadata$n_cells, 10))
  rhdf5::h5createDataset(file = h5_fp, dataset = "feature_idxs", dims = mtx_metadata$n_data_points, storage.mode = "integer", level = 0L, chunk = min(mtx_metadata$n_data_points - 1, 50))
  if (!mtx_metadata$is_logical) {
    rhdf5::h5createDataset(file = h5_fp, dataset = "data_csc", dims = mtx_metadata$n_data_points, storage.mode = "integer", level = 0L, chunk = min(mtx_metadata$n_data_points - 1, 50))
  }

  # Initialize CSR
  rhdf5::h5createDataset(file = h5_fp, dataset = "feature_ptr", dims = mtx_metadata$n_features + 1, storage.mode = "integer", level = 0L, chunk = min(mtx_metadata$n_features, 10))
  rhdf5::h5createDataset(file = h5_fp, dataset = "cell_idxs", dims = mtx_metadata$n_data_points, storage.mode = "integer", level = 0L, chunk = min(mtx_metadata$n_data_points - 1, 50))
  if (!mtx_metadata$is_logical) {
    rhdf5::h5createDataset(file = h5_fp, dataset = "data_csr", dims = mtx_metadata$n_data_points, storage.mode = "integer", level = 0L, chunk = min(mtx_metadata$n_data_points - 1, 50))
  }

  invisible(NULL)
}

#' Initialize accumulator
#'
#' @param terminal_symbol the terminal symbol
#' @param bag_of_variables the bag of variables
#'
#' @return An initialized vector of the correct type and length
#' @noRd
initialize_accumulator <- function(terminal_symbol, bag_of_variables) {
  info <- get_terminal_acc_and_args(terminal_symbol)
  info$acc_constructor(bag_of_variables[[info$acc_length]])
}


#' Get accumulator funct and arg list
#'
#' @param terminal_symbol a terminal symbol
#'
#' @return a list containing (i) the name of and (ii) the arguments to pass to the accumulator function.
#' @noRd
get_accumulator_funct_arg_list <- function(terminal_symbol) {
  info <- get_terminal_acc_and_args(terminal_symbol)
  list(acc_funct = info$acc_funct, acc_args = info$acc_args)
}


################################################################
# Core algorithm functions start
################################################################
#' Run mtx algo step 1 in: 1. chunk mode on one large .mtx file; 2. list of .mtx files mode
#'
#' Runs the first step of the .mtx algo.
#'
#' @param mtx_fp file path to .mtx file
#' @param initialize_accumulator initialize accumulator function
#' @param bag_of_variables the bag of variables
#' @param is_logical (boolean) is matrix logical
#' @param n_lines_per_chunk number of lines to process per chunk
#' @param n_rows_to_skip number of rows to skip in .mtx file
#'
#' @return list containing (i) n_features, and (ii) a list containing n_features an n_features vector for each chunk.
#' @noRd
run_core_algo_step_1 <- function(mtx_fp, initialize_accumulator, bag_of_variables, is_logical, n_lines_per_chunk, n_rows_to_skip, progress) {
  symbols <- symbols_enum()
  initializer <- function() initialize_accumulator(terminal_symbol = symbols$n_nonzero_feature, bag_of_variables = bag_of_variables)
  arguments <- arguments_enum()
  acc_init <- list(initializer(), list())

  # the function used to compute the accumulated row pointer
  get_accumulated_row_ptr <- function(x, pos, acc) {
    decrement_idxs(x$feature_idxs)
    n_nonzero_features_chunk <- initializer()
    inc_n_entries(n_nonzero_features_chunk, x$feature_idxs)
    sum_in_place(acc[[1]], n_nonzero_features_chunk)
    acc[[2]][[length(acc[[2]]) + 1]] <- n_nonzero_features_chunk
    return(acc)
  }

  if (length(mtx_fp) == 1) {
    return(run_core_algo_step_1_mtxchunked(mtx_fp, get_accumulated_row_ptr, bag_of_variables, is_logical, n_lines_per_chunk, n_rows_to_skip, progress, arguments, acc_init))
  } else {
    return(run_core_algo_step_1_mtxfilelist(mtx_fp, get_accumulated_row_ptr, bag_of_variables, is_logical, n_lines_per_chunk, n_rows_to_skip, progress, arguments, acc_init))
  }
}

#' Run mtx algo step 1 in list of .mtx files mode
#'
#' Runs the first step of the .mtx algo.
#'
#' @param mtx_fp file path to .mtx file
#' @param get_accumulated_row_ptr the function to compute accumulated row pointer
#' @param bag_of_variables the bag of variables
#' @param is_logical (boolean) is matrix logical
#' @param n_lines_per_chunk number of lines to process per chunk
#' @param n_rows_to_skip number of rows to skip in .mtx file
#' @param arguments arguments
#' @param acc_init empty accumulator, list of three elements: (i) a vector to store the n_nonzero_features accumulator, and (ii) a list to store the number of nonzero features in each chunk (iii) empty accumulator.
#'
#' @return list containing (i) n_features, and (ii) a list containing n_features an n_features vector for each chunk.
#' @noRd
run_core_algo_step_1_mtxfilelist <- function(mtx_fp, get_accumulated_row_ptr, bag_of_variables, is_logical, n_lines_per_chunk, n_rows_to_skip, progress, arguments, acc_init) {
  acc <- acc_init
  for (i in seq(1, length(mtx_fp))) {
    x <- readr::read_delim(file = mtx_fp[i],
                           delim = " ",
                           skip = n_rows_to_skip[i],
                           col_names = arguments$feature_idxs,
                           progress = progress,
                           col_types = if (is_logical) "i_" else "i__")
    acc <- get_accumulated_row_ptr(x, 0, acc)
  }
  return(acc)
}


#' Run mtx algo step 1 in: chunk mode on one large .mtx file
#'
#' Runs the first step of the .mtx algo.
#'
#' @param mtx_fp file path to .mtx file
#' @param get_accumulated_row_ptr the function to compute accumulated row pointer
#' @param bag_of_variables the bag of variables
#' @param is_logical (boolean) is matrix logical
#' @param n_lines_per_chunk number of lines to process per chunk
#' @param n_rows_to_skip number of rows to skip in .mtx file
#' @param arguments arguments
#' @param acc_init empty accumulator, list of three elements: (i) a vector to store the n_nonzero_features accumulator, and (ii) a list to store the number of nonzero features in each chunk (iii) empty accumulator.
#'
#' @return list containing (i) n_features, and (ii) a list containing n_features an n_features vector for each chunk.
#' @noRd
run_core_algo_step_1_mtxchunked <- function(mtx_fp, get_accumulated_row_ptr, bag_of_variables, is_logical, n_lines_per_chunk, n_rows_to_skip, progress,arguments, acc_init) {
  ret <- readr::read_delim_chunked(file = mtx_fp,
                            chunk_size = n_lines_per_chunk,
                            skip = n_rows_to_skip,
                            callback = readr::AccumulateCallback$new(get_accumulated_row_ptr, acc = acc_init),
                            delim = " ",
                            col_names = arguments$feature_idxs,
                            progress = progress,
                            col_types = if (is_logical) "i_" else "i__")
  return(ret)
}


#' Run subtask 2a
#'
#' Runs subtask a of part 2 of mtx algorithm: updates the accumulator of each terminal symbol.
#'
#' @param x a data.table (passed by ref)
#' @param bag_of_variables the bag_of_variables (also passed by ref)
#' @param acc the accumulator list
#' @param terminal_functs_args arguments to pass to the accumulator for a given terminal
#'
#' @return NULL
#' @noRd
run_subtask_2a <- function(x, bag_of_variables, acc, terminal_functs_args) {
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


#' Run subtask 2b
#'
#' Runs subtask 2b of algorithm: Writes chunk of gene_idxs and (if integer matrix) data to CSC matrix on-disk.
#'
#' NOTE: parts of this function should be rewritten in C++; in particular, we may need to handle pos manually, as r integers have limited size.
#'
#' @param x a data.table
#' @param pos the starting row in the data.table; uses 1-based indexing
#' @param h5_fp file path to the on-disk h5 file.
#' @param is_logical is the mtx logical
#'
#' @return NULL
#' @noRd
run_subtask_2b <- function(x, pos, h5_fp, is_logical) {
  # Write feature idxs
  write_data_h5(h5_fp, "feature_idxs", x$feature_idxs, pos - 1L)
  if (!is_logical) {
    # If integer matrix, write data too
    write_data_h5(h5_fp, "data_csc", x$umi_counts, pos - 1L)
  }
  return(invisible())
}


#' Run subtask 2c
#'
#' @param x a data.table
#' @param h5_fp file path to on-disk h5 file
#' @param is_logical (boolean) is the matrix logical?
#' @param row_ptr the (accumulated) row pointer
#' @param n_nonzero_features_per_chunk a list of the number of nonzero features in each chunk
#' @param chunk_no the current chunk number
#' @param n_features total number of features in matrix
#' @noRd
run_subtask_2c <- function(x, h5_fp, is_logical, row_ptr, n_nonzero_features_per_chunk, chunk_no, n_features) {
  arguments <- arguments_enum()
  data.table::setorderv(x, arguments$feature_idxs)
  n_nonzero_features_chunk <- n_nonzero_features_per_chunk[[chunk_no]]
  in_memory_row_ptr <- c(0L, cumsum(n_nonzero_features_chunk))
  if (!is_logical) {
    map_memory_to_disk(file_name_in = h5_fp,
                       m_cell_idxs = x$cell_idxs,
                       cell_idxs_name = "cell_idxs",
                       m_umi_counts = x$umi_counts,
                       umi_counts_name = "data_csr",
                       n_features = n_features,
                       m_row_ptr = in_memory_row_ptr,
                       f_row_ptr = row_ptr)
    } else {
      map_memory_to_disk_logical_matrix(file_name_in = h5_fp,
                                        m_cell_idxs = x$cell_idxs,
                                        cell_idxs_name = "cell_idxs",
                                        n_features = n_features,
                                        m_row_ptr = in_memory_row_ptr,
                                        f_row_ptr = row_ptr)
  }
  sum_in_place(row_ptr, n_nonzero_features_chunk)
}


#' Run mtx algo step 2
#'
#' This function runs step 2 of the core mtx algorithm. It (a) computes the terminal symbols, (b) writes to the CSC matrix, and (c) sorts the data by feature_idx, then writes to the CSR matrix.
#'
#' @param h5_fp full path to the h5 file on-disk
#' @param mtx_fp file path to mtx
#' @param is_logical is the mtx file logical
#' @param bag_of_variables the bag of variables containing the variables to pass to the accumulator functions
#' @param n_lines_per_chunk number of elements to load per chunk
#' @param n_rows_to_skip n rows to skip in mtx file
#' @param initial_accumulators list of starting accumulators
#' @param terminal_functs_args list of accumulator function names and arguments
#' @param row_ptr the starting row pointer
#' @param n_nonzero_features_per_chunk initial vector for n_nonzero_features_per_chunk
#'
#' @return a list containing the values of the terminals
#' @noRd
run_core_algo_step_2 <- function(h5_fp, mtx_fp, is_logical, bag_of_variables, n_lines_per_chunk, n_rows_to_skip, initial_accumulators, terminal_functs_args, row_ptr, n_nonzero_features_per_chunk, progress) {
  chunk_no <- 1L
  # Define closure to be called by readr::read_delim_chunked
  closure <- function(x, pos, acc) {
    # example chunk: x <- read.table(file = mtx_fp, header = FALSE, sep = " ", col.names = c("feature_idxs", "cell_idxs", if (is_logical) NULL else "umi_counts"), skip = n_rows_to_skip, colClasses = rep("integer", if (is_logical) 2 else 3), nrows = n_lines_per_chunk); pos <- 1
    data.table::setDT(x)
    decrement_idxs(x$feature_idxs)
    decrement_idxs(x$cell_idxs)

    # run subtask a
    run_subtask_2a(x, bag_of_variables, acc[[1]], terminal_functs_args)
    # run subtask b
    if (progress) cat("\nWriting CSC data.\n")
    run_subtask_2b(x, pos, h5_fp, is_logical)
    # run subtask c
    if (progress) cat("Writing CSR data.\n")
    run_subtask_2c(x, h5_fp, is_logical, acc[[2]], n_nonzero_features_per_chunk, chunk_no, bag_of_variables$n_features)
    # increment chunk_no in enclosing environment
    chunk_no <<- chunk_no + 1L
    return(acc)
  }
  arguments <- arguments_enum()
  # first element: terminal accumulator list; second element: row_ptr.
  acc_init <- list(initial_accumulators, row_ptr)

  if (length(mtx_fp) == 1) {
    return(run_core_algo_step_2_mtxchunked(mtx_fp, is_logical, n_lines_per_chunk, n_rows_to_skip, progress, closure, acc_init, arguments))
  } else {
    return(run_core_algo_step_2_mtxfilelist(mtx_fp, is_logical, n_lines_per_chunk, n_rows_to_skip, progress, closure, acc_init, arguments, bag_of_variables))
  }
}

#' Run mtx algo step 2 in: chunk mode on one large .mtx file
#'
#' @param mtx_fp file path to mtx
#' @param is_logical is the mtx file logical
#' @param n_lines_per_chunk number of elements to load per chunk
#' @param n_rows_to_skip n rows to skip in mtx file
#' @param progress progress
#' @param closure closure function to calculate the terminals and write to the h5 file
#' @param acc_init empty accumulator
#' @param arguments arguments
#'
#' @return list containing the values of the terminals
#' @noRd
run_core_algo_step_2_mtxchunked <- function(mtx_fp, is_logical, n_lines_per_chunk, n_rows_to_skip, progress, closure, acc_init, arguments) {
  terminals <- readr::read_delim_chunked(file = mtx_fp,
                                         chunk_size = n_lines_per_chunk,
                                         skip = n_rows_to_skip,
                                         callback = readr::AccumulateCallback$new(closure, acc = acc_init),
                                         delim = " ",
                                         col_names = c(arguments$feature_idxs, arguments$cell_idxs, if (is_logical) NULL else arguments$umi_counts),
                                         progress = progress,
                                         col_types = if (is_logical) "ii" else "iii")
  return(terminals)
}

#' Run mtx algo step 2 in list of .mtx files mode
#'
#' @param mtx_fp file path to mtx
#' @param is_logical is the mtx file logical
#' @param n_lines_per_chunk number of elements to load per chunk
#' @param n_rows_to_skip n rows to skip in mtx file
#' @param progress progress
#' @param closure closure function to calculate the terminals and write to the h5 file
#' @param acc_init empty accumulator
#' @param arguments arguments
#' @param bag_of_variables the bag of variables containing the variables to pass to the accumulator functions
#'
#' @return list containing the values of the terminals
#' @noRd
run_core_algo_step_2_mtxfilelist <- function(mtx_fp, is_logical, n_lines_per_chunk, n_rows_to_skip, progress, closure, acc_init, arguments, bag_of_variables) {
  acc <- acc_init
  cell_idx <- as.integer(c(0, cumsum(bag_of_variables$n_cells_in_files)))
  pos <- 1L
  for (i in seq(1, length(mtx_fp))) {
    x <- readr::read_delim(file = mtx_fp[i],
                           delim = " ",
                           skip = n_rows_to_skip[i],
                           col_names = c(arguments$feature_idxs, arguments$cell_idxs, if (is_logical) NULL else arguments$umi_counts),
                           progress = progress,
                           col_types = if (is_logical) "ii" else "iii")
    #increase col_idx
    x$cell_idxs = x$cell_idxs + cell_idx[i]
    acc <- closure(x, pos, acc)
    pos <- pos + nrow(x)
  }
  return(acc)
}

#' Run core mtx algo
#'
#' Runs to core algorithm for mtx files. There are two steps:
#' (i) compute the row pointer
#' (ii) write CSC matrix, write CSR matrix, compute covariate matrices
#'
#' @param h5_fp file path to h5 file
#' @param mtx_fp file path to .mtx file
#' @param is_logical (logical) is expression matrix logical?
#' @param covariates a list of two elements: the feature covariates and the cell covariates
#' @param bag_of_variables environment containing named variables
#' @param n_lines_per_chunk number of elements to load per chunk
#' @param n_rows_to_skip number of rows to skip in .mtx file
#' @noRd
run_core_algo <- function(h5_fp, mtx_fp, is_logical, covariates, bag_of_variables, n_lines_per_chunk, n_rows_to_skip, progress) {
  grammar <- initialize_grammar()
  symbols <- symbols_enum()

  # Run step 1 of core algorithm
  row_ptrs <- run_core_algo_step_1(mtx_fp, initialize_accumulator, bag_of_variables,
                                  is_logical, n_lines_per_chunk, n_rows_to_skip, progress)

  # Compute row pointer and write to disk.
  row_ptr <- c(0L, cumsum(row_ptrs[[1]]))
  write_data_h5(file_name_in = h5_fp, dataset_name_in = "feature_ptr", buffer = row_ptr, start_pos = 0)

  # Prepare step 2 of core algo; determine which terminal symbols to compute
  terminal_symbols <- lapply(unlist(covariates),
                             get_terminals_for_covariate, grammar = grammar) %>% unlist() %>% unique()

  # Exclude from this vector n_nonzero_features, which we already computed in step 1.
  terminal_symbols <- terminal_symbols[!(terminal_symbols == symbols$n_nonzero_feature)]

  # Obtain the starting accumulator, as well as the accumulator function name and args, for each terminal
  initial_accumulators <- lapply(terminal_symbols, initialize_accumulator, bag_of_variables = bag_of_variables)
  terminal_functs_args <- lapply(terminal_symbols, get_accumulator_funct_arg_list)

  # Run step 2 of core algorithm
  terminal_values_and_row_ptr <- run_core_algo_step_2(h5_fp, mtx_fp, is_logical, bag_of_variables,
                                                     n_lines_per_chunk, n_rows_to_skip, initial_accumulators,
                                                     terminal_functs_args, row_ptr, row_ptrs[[2]], progress)
  # Compute and write column pointer to CSC matrix
  terminal_values <- terminal_values_and_row_ptr[[1]]
  n_nonzero_cell <- terminal_values[[which(terminal_symbols == symbols$n_nonzero_cell)]]
  cell_ptr <- c(0, cumsum(n_nonzero_cell))
  rhdf5::h5write(cell_ptr, h5_fp, "cell_ptr")

  # compute the covariate matrices
  for (i in 1:length(terminal_symbols)) grammar[[terminal_symbols[i]]]$value <- terminal_values[[i]]
  grammar[[symbols$n_nonzero_feature]]$value <- row_ptrs[[1]]
  cov_mats <- lapply(covariates, function(covariate_vect) {
    out <- lapply(covariate_vect, evaluate_grammar, grammar) %>% as.data.frame()
    colnames(out) <- gsub('_(feature|cell)', "", covariate_vect)
    return(out)
  })

  # return covariate matrices
  return(cov_mats)
}
