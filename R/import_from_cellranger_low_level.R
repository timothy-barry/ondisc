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


initialize_cellwise_covariates <- function(modality_names, n_cells, compute_cell_cycle = FALSE) {
  possible_covariates <- c("n_umis", "n_nonzero", "p_mito", "feature_w_max_expression", "frac_umis_max_feature", "s_score", "g2m_score", "phase")
  covariates_to_compute <- list("Gene Expression" = c("n_umis", "n_nonzero", "p_mito"),
                                "CRISPR Guide Capture" = c("n_umis", "n_nonzero", "feature_w_max_expression", "frac_umis_max_feature"),
                                "Antibody Capture" = c("n_umis", "n_nonzero"))

  # Add cell cycle covariates to Gene Expression if requested
  if (compute_cell_cycle && "Gene Expression" %in% modality_names) {
    covariates_to_compute[["Gene Expression"]] <- c(covariates_to_compute[["Gene Expression"]], "s_score", "g2m_score")
  }
  out <- lapply(modality_names, function(modality_name) {
    # evaluate the boolean vector for the modality
    curr_covariates <- covariates_to_compute[[modality_name]]
    bool_vect <- vapply(possible_covariates, function(possible_covariate) {
      possible_covariate %in% curr_covariates
    }, FUN.VALUE = logical(1))
    covariate_list <- list(
      n_umis = integer(length = if (bool_vect[["n_umis"]]) n_cells else 0L),
      n_nonzero = integer(length = if (bool_vect[["n_nonzero"]]) n_cells else 0L),
      p_mito = numeric(length = if (bool_vect[["p_mito"]]) n_cells else 0L),
      feature_w_max_expression = integer(length = if (bool_vect[["feature_w_max_expression"]]) n_cells else 0L),
      frac_umis_max_feature = numeric(length = if (bool_vect[["frac_umis_max_feature"]]) n_cells else 0L),
      s_score = numeric(length = if (bool_vect[["s_score"]]) n_cells else 0L),
      g2m_score = numeric(length = if (bool_vect[["g2m_score"]]) n_cells else 0L),
      phase = if (bool_vect[["phase"]]) factor(rep(NA, n_cells), levels = c("G1", "S", "G2M", "Undecided")) else factor()
    )
    list(bool_vect = bool_vect,
         covariate_list = covariate_list)
  }) |> stats::setNames(modality_names)
  return(out)
}


prepare_output_covariate_dt <- function(cellwise_covariates, new_modality_names, n_cells_per_batch, modality_feature_ids) {
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
                            modality_feature_ids, n_cells, integer_id, chunk_size, compression_level) {
  row_ptr_list <- lapply(seq_along(modality_names), function(k) {
    file_name_in <- file_names_in[k]
    create_odm(file_name_in = file_name_in,
               n_nonzero_features = n_nonzero_features_vector_list[[k]],
               feature_ids = modality_feature_ids[[k]],
               n_cells = n_cells,
               integer_id = integer_id,
               chunk_size = chunk_size,
               compression_level = compression_level)
    read_row_ptr(file_name_in, length(n_nonzero_features_vector_list[[k]]))
  }) |> stats::setNames(modality_names)
  return(row_ptr_list)
}


process_input_files_round_1 <- function(matrix_fps, modality_names, modality_start_idx_features,
                                        feature_idx_to_vector_idx_map, n_features_per_modality,
                                        compute_cell_cycle = FALSE, cc_scale_factor = 10000) {
  # 0. initialize the outputs:
  # a. modality_start_idx_mtx_list
  # b. n_nonzero_features_vector_list
  # c. n_cells
  # d. n_cells_per_batch
  # e. collapse_grna_counts
  modality_start_idx_mtx_list <- vector(mode = "list", length = length(matrix_fps))
  n_nonzero_features_vector_list <- lapply(seq_along(modality_names), function(k) {
    integer(length = n_features_per_modality[k])
  }) |> stats::setNames(modality_names)
  n_cells_per_batch <- integer(length = length(matrix_fps))
  collapse_grna_counts <- !is.null(feature_idx_to_vector_idx_map)

  # Initialize gene_norm_sum if cell cycle scoring is requested
  gene_norm_sum <- NULL
  if (compute_cell_cycle && "Gene Expression" %in% modality_names) {
    ge_idx <- which(modality_names == "Gene Expression")
    n_ge_features <- n_features_per_modality[ge_idx]
    gene_norm_sum <- rep(0.0, n_ge_features)
  }

  cat("Round 1/2 processing of the input files.\n")
  for (i in seq_along(matrix_fps)) {
    cat(paste0("\tProcessing file ", i , " of ", length(matrix_fps), ".\n"))
    # 1. set the matrix fp and obtain the metadata
    matrix_fp <- matrix_fps[i]
    mtx_metadata <- get_mtx_metadata(matrix_fp)
    n_cells_per_batch[i] <- mtx_metadata$n_cells

    # 2. load the matrix data - determine columns based on needed functionality
    need_cell_idx <- collapse_grna_counts
    need_umi_counts <- compute_cell_cycle && "Gene Expression" %in% modality_names

    if (need_umi_counts) {
      # Read all three columns when cell cycle scoring is enabled
      my_col_names <- c("feature_idx", "j", "x")
      my_col_classes <- c("integer", "integer", "integer")
    } else if (need_cell_idx) {
      # Read feature and cell indices when gRNA collapsing is needed
      my_col_names <- c("feature_idx", "j")
      my_col_classes <- c("integer", "integer", "NULL")
    } else {
      # Read only feature indices for basic processing
      my_col_names <- "feature_idx"
      my_col_classes <- c("integer", "NULL", "NULL")
    }

    dt <- data.table::fread(file = matrix_fp,
                            skip = mtx_metadata$n_to_skip,
                            col.names = my_col_names,
                            colClasses = my_col_classes,
                            showProgress = FALSE, nThread = 1)
    decrement_vector(dt$feature_idx) # decrement feature idx
    data.table::setkey(dt, feature_idx) # radix sort on feature_idx

    # 3. determine the start idx of the different modalities; convert into a list
    modality_start_mtx <- c(vapply(seq_along(modality_names), function(k) {
      modality_start_idx <- modality_start_idx_features[[k]]
      dt[list(modality_start_idx), which = TRUE, mult = "first", roll = -Inf] - 1L # binary search for start idxs
    }, FUN.VALUE = integer(1)), nrow(dt))
    modality_start_mtx <- lapply(seq_along(modality_names), function(k) {
      c(modality_start_mtx[k], modality_start_mtx[k + 1L] - 1L)
    }) |> stats::setNames(modality_names)
    modality_start_idx_mtx_list[[i]] <- modality_start_mtx

    # 4. (possibly) collapse grna counts
    if (collapse_grna_counts) {
      modality_start_mtx <- collapse_grna_counts(dt = dt,
                                                 feature_idx_to_vector_idx_map = feature_idx_to_vector_idx_map,
                                                 modality_start_mtx = modality_start_mtx,
                                                 round_1 = TRUE)
    }

    # 5. (optionally) collect gene expression statistics for cell cycle scoring
    if (compute_cell_cycle && "Gene Expression" %in% modality_names) {
      # Find Gene Expression modality range
      ge_idx <- which(modality_names == "Gene Expression")
      ge_start_idx <- modality_start_mtx[[ge_idx]][1]
      ge_end_idx <- modality_start_mtx[[ge_idx]][2]

      if (ge_end_idx >= ge_start_idx) {
        # First, compute UMI totals per cell for this batch
        current_n_cells <- n_cells_per_batch[i]
        batch_n_umis <- integer(current_n_cells)

        # Accumulate UMI totals for Gene Expression genes only
        # Note: ge_start_idx and ge_end_idx are 0-based indices into the data table
        for (row_idx in (ge_start_idx + 1):(ge_end_idx + 1)) {  # Convert to 1-based for R indexing
          cell_idx <- dt$j[row_idx] + 1  # Convert 0-based cell index to 1-based
          batch_n_umis[cell_idx] <- batch_n_umis[cell_idx] + dt$x[row_idx]
        }

        # Compute normalized gene expression statistics
        compute_normalized_gene_expression_stats(
          gene_norm_sum = gene_norm_sum,
          feature_idx = dt$feature_idx,
          j = dt$j,
          x = dt$x,
          n_umis = batch_n_umis,
          start_idx = ge_start_idx,
          end_idx = ge_end_idx,
          feature_offset = modality_start_idx_features[[ge_idx]],
          scale_factor = cc_scale_factor
        )
      }
    }

    # 6. update n_nonzero_features_vector for each modality
    for (k in seq_along(modality_names)) {
      start_idx <- modality_start_mtx[[k]][1]
      end_idx <- modality_start_mtx[[k]][2]
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
              n_cells_per_batch = n_cells_per_batch,
              gene_norm_sum = gene_norm_sum))
}


process_input_files_round_2 <- function(matrix_fps, file_names_in, modality_names, modality_start_idx_features,
                                        row_ptr_list, modality_start_idx_mtx_list, mt_feature_idxs,
                                        cellwise_covariates, feature_idx_to_vector_idx_map,
                                        compute_cell_cycle = FALSE, cc_gene_indices = NULL, cc_control_genes = NULL, cc_scale_factor = 10000) {
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

    # 3. collapse grna counts
    modality_start_mtx <- modality_start_idx_mtx_list[[i]]
    if (!is.null(feature_idx_to_vector_idx_map)) {
      modality_start_mtx <- collapse_grna_counts(dt = dt,
                                                 feature_idx_to_vector_idx_map = feature_idx_to_vector_idx_map,
                                                 modality_start_mtx = modality_start_mtx,
                                                 round_1 = FALSE)
    }

    # 4. compute the cell-wise covariates for each modality
    cat("Computing cellwise covariates. ")
    for (k in seq_along(modality_names)) {
      start_idx <- modality_start_mtx[[k]][1]
      end_idx <- modality_start_mtx[[k]][2]
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

    # 4.5. Compute cell cycle scores if requested (only for Gene Expression modality)
    if (compute_cell_cycle && "Gene Expression" %in% modality_names) {
      cat("Computing cell cycle scores. ")
      ge_idx <- which(modality_names == "Gene Expression")
      start_idx <- modality_start_mtx[[ge_idx]][1]
      end_idx <- modality_start_mtx[[ge_idx]][2]
      feature_offset <- modality_start_idx_features[[ge_idx]]

      compute_cell_cycle_scores(
        s_scores = cellwise_covariates[[ge_idx]]$covariate_list$s_score,
        g2m_scores = cellwise_covariates[[ge_idx]]$covariate_list$g2m_score,
        s_gene_indices = cc_gene_indices$s_gene_indices,
        g2m_gene_indices = cc_gene_indices$g2m_gene_indices,
        s_control_indices = cc_control_genes$s_control_indices,
        g2m_control_indices = cc_control_genes$g2m_control_indices,
        feature_idx = dt$feature_idx,
        j = dt$j,
        x = dt$x,
        n_umis = cellwise_covariates[[ge_idx]]$covariate_list$n_umis,
        start_idx = start_idx,
        end_idx = end_idx,
        feature_offset = feature_offset,
        cell_offset = n_cum_cells,
        n_cells = n_cells,
        scale_factor = cc_scale_factor
      )
      cat("Done.\n")
    }

    # 5. Write the data to disk in CSR format
    cat("Writing to disk.\n")
    for (k in seq_along(modality_names)) {
      start_idx <- modality_start_mtx[[k]][1]
      end_idx <- modality_start_mtx[[k]][2]
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

    # 6. increment n_cum_cells
    n_cum_cells <- n_cum_cells + n_cells
    rm(dt); gc() |> invisible()
  }
  return(NULL)
}
