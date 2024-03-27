my_seed <- as.integer(Sys.time() |> format("%H%M%S"))
print(paste0("seed: ", my_seed))
set.seed(my_seed)

test_that("import data from cellranger v2", {
  ########################
  # define test parameters
  ########################
  n_trials <- 5L
  n_rows_range <- c(10L, 10000L)
  n_cols_range <- c(10L, 10000L)
  max_mito_genes <- 0.1
  max_n_batches <- 5
  # trial_w_vector <- sample(x = seq(1, n_trials), size = 1L, replace = FALSE)
  trial_w_vector <- n_trials + 1L # do not use vector for now

  #################
  # create the data
  #################
  test_data_list <- lapply(seq(1L, n_trials), FUN = function(i) {
    print(paste0("Generating example dataset ", i, "."))
    modalities <- sample(x = c("gene", "grna", "protein"),
                         size = sample(seq(1, 3), 1), replace = FALSE)
    if (i == trial_w_vector) modalities <- c("grna", modalities) |> unique()
    n_features <- sample(x = seq(n_rows_range[1], n_rows_range[2]),
                         size = length(modalities), replace = FALSE)
    n_cells <- sample(x = seq(n_cols_range[1], n_cols_range[2]), size = 1L)
    n_batch <- sample(x = seq(1, max_n_batches), size = 1)
    curr_base_directory <- create_new_directory()
    curr_data <- write_example_cellranger_dataset(n_features = n_features,
                                                  n_cells = n_cells,
                                                  n_batch = n_batch,
                                                  modalities = modalities,
                                                  dir_to_write = curr_base_directory,
                                                  p_zero = min(runif(1), 0.9),
                                                  p_set_col_zero = min(runif(1), 0.9),
                                                  p_set_row_zero = min(runif(1), 0.9))
    # check to make sure some entries are nonzero; if so, return
    if (any(sapply(curr_data$matrix_list, function(mat) length(mat@x)) == 0L)) {
      return()
    }
    # construct vector data frame if trials_w_vector
    if (i == trial_w_vector) {
      grna_matrix <- curr_data$matrix_list$grna
      n_vectors <- ceiling(nrow(grna_matrix)/8)
      sample_probs <- runif(n_vectors)
      sample_probs <- sample_probs/sum(sample_probs)
      vector_idxs <- sample(x = seq(1L, n_vectors), size = nrow(grna_matrix),
                            replace = TRUE, prob = sample_probs)
      vector_ids <- paste0("vector-", vector_idxs)
      grna_target_data_frame <- data.frame(grna_id = rownames(grna_matrix),
                                           vector_id = vector_ids)
    } else {
      grna_target_data_frame <- NULL
    }

    curr_data$grna_target_data_frame <- grna_target_data_frame
    curr_data$base_directory <- curr_base_directory
    return(curr_data)
  })

  ###############
  # run the tests
  ###############
  for (i in seq(1L, n_trials)) {
    print(paste0("Testing import from cellranger for dataset ", i))
    directories_to_load <- list.files(test_data_list[[i]]$base_directory,
                                      full.names = TRUE, pattern = "batch_*")
    # load the data from cellranger format
    odms <- create_odm_from_cellranger(directories_to_load = directories_to_load,
                                       directory_to_write = test_data_list[[i]]$base_directory,
                                       grna_target_data_frame = test_data_list[[i]]$grna_target_data_frame,
                                       chunk_size = 5L)

    # loop over the moalities, testing dimension and load functionality
    modalities <- names(test_data_list[[i]]$matrix_list)
    for (modality in modalities) {
      odm <- odms[[modality]]
      mem_matrix <- test_data_list[[i]]$matrix_list[[modality]]

      # if grna modality, perform additional processing
      if (modality == "grna") {
        if (i == trial_w_vector) {
          grna_target_data_frame <- test_data_list[[i]]$grna_target_data_frame
          vector_ids <- unique(grna_target_data_frame$vector_id)
          mem_matrix <- sapply(vector_ids, function(curr_vector_id) {
            curr_grna_ids <- grna_target_data_frame |>
              dplyr::filter(vector_id == curr_vector_id) |>
              dplyr::pull(grna_id)
            Matrix::colSums(mem_matrix[curr_grna_ids,,drop=FALSE])
          }) |> t()
          test_data_list[[i]]$grna_matrix <- mem_matrix
        }
      }
      # 1. check dimension
      expect_equal(dim(mem_matrix), odm@dimension)
      # 2. check feature ids
      expect_equal(rownames(mem_matrix), rownames(odm))
      # 3. check index into randomly selected rows by integer
      n_row <- nrow(mem_matrix)
      sample_idxs <- c(1L, sample(x = seq(1L, n_row), size = min(30, n_row), replace = FALSE), n_row)
      for (sample_idx in sample_idxs) {
        expect_equal(odm[sample_idx,], as.integer(mem_matrix[sample_idx,]))
      }
      # 4. check index into randomly selected rows by feature
      sample_features <- sample(x = rownames(mem_matrix), size = min(30, n_row), replace = FALSE)
      for (sample_feature in sample_features) {
        expect_equal(odm[sample_feature,], as.integer(mem_matrix[sample_feature,]))
      }
      # 5. check the covariates
      cellwise_covariates_disk <- odms$cellwise_covariates
      covariates_to_compute <- list("gene" = c("n_umis", "n_nonzero", "p_mito"),
                                    "grna" = c("n_umis", "n_nonzero", "feature_w_max_expression", "frac_umis_max_feature"),
                                    "protein" = c("n_umis", "n_nonzero"))
      covariates_to_check <- covariates_to_compute[[modality]]
      for (curr_covariate in covariates_to_check) {
        if (curr_covariate == "n_umis") {
          n_umis_mem <- Matrix::colSums(mem_matrix)
          n_umis_disk <- cellwise_covariates_disk[[paste0(modality, "_n_umis")]]
          expect_equal(n_umis_mem, n_umis_disk)
        } else if (curr_covariate == "n_nonzero") {
          n_nonzero_mem <- Matrix::colSums(mem_matrix >= 1)
          n_nonzero_disk <- cellwise_covariates_disk[[paste0(modality, "_n_nonzero")]]
          expect_equal(n_umis_mem, n_umis_disk)
        } else if (curr_covariate == "p_mito") {
          gene_names <- test_data_list[[i]]$gene_names
          n_umis_mito_mem <- Matrix::colSums(mem_matrix[grep(pattern = "^MT-", x = gene_names),])
          p_mito_mem <- ifelse(n_umis_mem == 0, 0, n_umis_mito_mem/n_umis_mem)
          p_mito_disk <- cellwise_covariates_disk$gene_p_mito
          expect_equal(p_mito_disk, p_mito_mem)
        } else if (curr_covariate == "feature_w_max_expression") {
          feature_idx_w_max_expression <- apply(X = mem_matrix, MARGIN = 2, FUN = which.max)
          feature_w_max_expression_mem <- rownames(mem_matrix)[feature_idx_w_max_expression]
          feature_w_max_expression_disk <- cellwise_covariates_disk$grna_feature_w_max_expression
          expect_equal(feature_w_max_expression_mem, feature_w_max_expression_disk)
        } else if (curr_covariate == "frac_umis_max_feature") {
          max_exp <- apply(X = mem_matrix, MARGIN = 2, FUN = max)
          frac_umis_max_feature <- ifelse(max_exp == 0, 1, max_exp / n_umis_mem)
          frac_umis_max_feature_disk <- cellwise_covariates_disk$grna_frac_umis_max_feature
          expect_equal(frac_umis_max_feature, frac_umis_max_feature_disk)
        }
      }
    }
  }
})
