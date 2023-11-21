get_mtx_metadata <- function(mtx_file, n_cells_col_id = 2L) {
  con <- file(mtx_file, "r")
  n_to_skip <- 0L
  while (TRUE) {
    n_to_skip <- n_to_skip + 1L
    line <- readLines(con, n = 1)
    if (substr(line, 0, 1) != "%") break()
  }
  close(con)
  metrics <- strsplit(line, split = " ", fixed = TRUE)[[1]] |> as.integer()
  out <- list(n_cells = metrics[2], n_to_skip = n_to_skip)
  return(out)
}


initialize_cellwise_covariates <- function(modality_names, n_cells) {
  possible_covariates <- c("n_umis", "n_nonzero", "p_mito", "feature_w_max_expression", "frac_umis_max_feature")
  covariates_to_compute <- list("Gene Expression" = c("n_umis", "n_nonzero", "p_mito"),
                                "CRISPR Guide Capture" = c("n_umis", "n_nonzero", "feature_w_max_expression", "frac_umis_max_feature"))
  out <- lapply(modality_names, function(modality_name) {
    # evaluate the boolean vector for the modality
    bool_vect <- if (modality_name %in% names(covariates_to_compute)) {
      curr_covariates <- covariates_to_compute[[modality_name]]
      sapply(possible_covariates, function(possible_covariate) possible_covariate %in% curr_covariates)
    } else {
      rep(TRUE, length(possible_covariates)) |> stats::setNames(possible_covariates)
    }
    covariate_list <- list(
      n_umis = integer(length = if (bool_vect[["n_umis"]]) n_cells else 0L),
      n_nonzero = integer(length = if (bool_vect[["n_nonzero"]]) n_cells else 0L),
      p_mito = numeric(length = if (bool_vect[["p_mito"]]) n_cells else 0L),
      feature_w_max_expression = integer(length = if (bool_vect[["feature_w_max_expression"]]) n_cells else 0L),
      frac_umis_max_feature = numeric(length = if (bool_vect[["frac_umis_max_feature"]]) n_cells else 0L)
    )
    list(bool_vect = bool_vect,
         covariate_list = covariate_list)
  }) |> stats::setNames(modality_names)
  return(out)
}


preprare_output_covariate_dt <- function(cellwise_covariates, new_modality_names, n_cells_per_batch, modality_feature_ids) {
  names(cellwise_covariates) <- new_modality_names
  l <- list()
  for (i in seq_along(new_modality_names)) {
    modality_name <- new_modality_names[i]
    computed_feats <- cellwise_covariates[[modality_name]]$bool_vect
    for (feat in names(computed_feats)) {
      if (computed_feats[[feat]]) {
        v <- cellwise_covariates[[modality_name]]$covariate_list[[feat]]
        if (feat == "feature_w_max_expression") {
          feature_ids <- modality_feature_ids[[i]]
          v <- feature_ids[v + 1L]
        }
        l[[paste0(modality_name, "_", feat)]] <- v
      }
    }
  }
  batch <- paste0("batch_", rep(x = seq_along(n_cells_per_batch), times = n_cells_per_batch)) |> factor()
  l$batch <- batch
  data.table::setDT(l)
}


initialize_odms <- function(modality_names, file_names_in, n_nonzero_features_vector_list,
                            modality_feature_ids, n_cells, chunk_size, compression_level) {
  row_ptr_list <- lapply(seq_along(modality_names), function(k) {
    file_name_in <- file_names_in[k]
    create_odm(file_name_in = file_name_in,
               n_nonzero_features = n_nonzero_features_vector_list[[k]],
               feature_ids = modality_feature_ids[[k]],
               n_cells = n_cells,
               chunk_size = chunk_size,
               compression_level = compression_level)
    read_row_ptr(file_name_in, length(n_nonzero_features_vector_list[[k]]))
  }) |> stats::setNames(modality_names)
  return(row_ptr_list)
}


process_input_files_round_1 <- function(matrix_fps, modality_names, modality_start_idx_features) {
  # 0. initialize the outputs:
  # a. modality_start_idx_mtx_list
  # b. n_nonzero_features_vector_list
  # c. n_cells
  # d. n_cells_per_batch
  n_features_per_modality <- diff(modality_start_idx_features)
  modality_start_idx_mtx_list <- vector(mode = "list", length = length(matrix_fps))
  n_nonzero_features_vector_list <- lapply(seq_along(modality_names), function(k) {
    integer(length = n_features_per_modality[k])
  }) |> stats::setNames(modality_names)
  n_cells_per_batch <- integer(length = length(matrix_fps))

  cat("Round 1/2 processing of the input files.\n")
  for (i in seq_along(matrix_fps)) {
    cat(paste0("\tProcessing file ", i , " of ", length(matrix_fps), ".\n"))
    # 1. set the matrix fp and obtain the metadata
    matrix_fp <- matrix_fps[i]
    mtx_metadata <- get_mtx_metadata(matrix_fp)
    n_cells_per_batch[i] <- mtx_metadata$n_cells

    # 2. load the feature idxs, decrement, and sort
    dt <- data.table::fread(file = matrix_fp,
                            skip = mtx_metadata$n_to_skip,
                            col.names = c("feature_idx"),
                            colClasses = c("integer", "NULL", "NULL"),
                            showProgress = FALSE, nThread = 1)
    decrement_vector(dt$feature_idx) # decrement
    data.table::setkey(dt, feature_idx) # radix sort on feature_idx

    # 3. determine the start idx of the different modalities
    modality_start_mtx <- sapply(seq_along(modality_names), function(k) {
      modality_start_idx <- modality_start_idx_features[[k]]
      dt[list(modality_start_idx), which = TRUE, mult = "first", roll = -Inf] - 1L # binary search for start idxs
    })
    modality_start_idx_mtx_list[[i]] <- modality_start_mtx <- c(modality_start_mtx, nrow(dt))

    # 4. update n_nonzero_features_vector for each modality
    for (k in seq_along(modality_names)) {
      start_idx <- modality_start_mtx[k]
      end_idx <- modality_start_mtx[k + 1]
      update_n_features_vector(feature_idx = dt$feature_idx,
                               n_nonzero_features_vector = n_nonzero_features_vector_list[[k]],
                               offset = modality_start_idx_features[[k]],
                               start_idx = start_idx,
                               end_idx = end_idx)
    }
    rm(dt); gc() |> invisible()
  }

  return(list(modality_start_idx_mtx_list = modality_start_idx_mtx_list,
              n_nonzero_features_vector_list = n_nonzero_features_vector_list,
              n_cells = sum(n_cells_per_batch),
              n_cells_per_batch = n_cells_per_batch))
}


process_input_files_round_2 <- function(matrix_fps, file_names_in, modality_names, modality_start_idx_features,
                                        row_ptr_list, modality_start_idx_mtx_list, mt_feature_idxs,
                                        cellwise_covariates) {
  cat("Round 2/2 processing of the input files.\n")
  n_cum_cells <- 0L
  n_features <- diff(modality_start_idx_features) |> stats::setNames(modality_names)
  for (i in seq_along(matrix_fps)) {
    cat(paste0("\tProcessing file ", i , " of ", length(matrix_fps), ". "))
    # 1. set the matrix fp and obtain the metadata
    matrix_fp <- matrix_fps[i]
    mtx_metadata <- get_mtx_metadata(matrix_fp)
    n_cells <- mtx_metadata$n_cells
    # 2. load the feature idxs, decrement, and sort
    dt <- data.table::fread(file = matrix_fp, skip = mtx_metadata$n_to_skip,
                            col.names = c("feature_idx", "j", "x"),
                            colClasses = c("integer", "integer", "integer"),
                            showProgress = FALSE, nThread = 1)
    decrement_vector(dt$feature_idx) # decrement feature_idx
    add_value_to_vector(dt$j, n_cum_cells - 1L)
    data.table::setorderv(dt, cols = c("feature_idx", "j")) # radix sort on feature_idx, cell_idx

    # 3. compute the cell-wise covariates for each modality
    cat("Computing cellwise covariates. ")
    for (k in seq_along(modality_names)) {
      start_idx <- modality_start_idx_mtx_list[[i]][k]
      end_idx <- modality_start_idx_mtx_list[[i]][k + 1]
      feature_offset <- modality_start_idx_features[[k]]
      compute_cellwise_covariates(n_umis = cellwise_covariates[[k]]$covariate_list$n_umis,
                                  n_nonzero = cellwise_covariates[[k]]$covariate_list$n_nonzero,
                                  p_mito = cellwise_covariates[[k]]$covariate_list$p_mito,
                                  feature_w_max_expression = cellwise_covariates[[k]]$covariate_list$feature_w_max_expression,
                                  frac_umis_max_feature = cellwise_covariates[[k]]$covariate_list$frac_umis_max_feature,
                                  feature_idx = dt$feature_idx,
                                  j = dt$j,
                                  x = dt$x,
                                  mt_feature_idxs = mt_feature_idxs[[k]],
                                  start_idx = start_idx,
                                  end_idx = end_idx,
                                  feature_offset = feature_offset,
                                  cell_offset = n_cum_cells,
                                  n_cells = n_cells,
                                  compute_n_umis = cellwise_covariates[[k]]$bool_vect[["n_umis"]],
                                  compute_n_nonzero = cellwise_covariates[[k]]$bool_vect[["n_nonzero"]],
                                  compute_p_mito = cellwise_covariates[[k]]$bool_vect[["p_mito"]],
                                  compute_feature_w_max_expression = cellwise_covariates[[k]]$bool_vect[["feature_w_max_expression"]],
                                  compute_frac_umis_max_feature = cellwise_covariates[[k]]$bool_vect[["frac_umis_max_feature"]])
    }

    # 4. Write the data to disk in CSR format
    cat("Writing to disk.\n")
    for (k in seq_along(modality_names)) {
      start_idx <- modality_start_idx_mtx_list[[i]][k]
      end_idx <- modality_start_idx_mtx_list[[i]][k + 1]
      feature_offset <- modality_start_idx_features[[k]]
      file_name_in <- file_names_in[k]
      write_to_csr(file_name_in = file_name_in,
                   start_idx = start_idx,
                   end_idx = end_idx,
                   feature_offset = feature_offset,
                   cell_offset = n_cum_cells,
                   n_features = n_features[k],
                   f_row_ptr = row_ptr_list[[k]],
                   feature_idx = dt$feature_idx,
                   m_j = dt$j,
                   m_x = dt$x)
    }

    # increment n_cum_cells
    n_cum_cells <- n_cum_cells + n_cells
    rm(dt); gc() |> invisible()
  }
  return(NULL)
}
