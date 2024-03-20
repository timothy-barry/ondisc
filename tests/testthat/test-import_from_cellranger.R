test_that("import data from cellranger", {
  ########################
  # define test parameters
  ########################
  n_trials <- 5L
  n_trials_with_vector_id <- 2L
  n_rows_range <- c(10L, 10000L)
  n_cols_range <- c(10L, 10000L)
  max_mito_genes <- 0.1
  max_n_batches <- 5
  trials_w_vector <- sample(x = seq(1, n_trials), size = n_trials_with_vector_id, replace = FALSE)

  #################
  # create the data
  #################
  test_data_list <- lapply(seq(1L, n_trials), FUN = function(i) {
    print(paste0("Generating example dataset ", i, "."))
    n_rows <- sample(x = seq(n_rows_range[1], n_rows_range[2]), size = 2L, replace = FALSE)
    n_col <- sample(x = seq(n_cols_range[1], n_cols_range[2]), size = 1L)
    gene_matrix <- create_random_matrix(n_row = n_rows[1],
                                        n_col = n_col,
                                        p_zero = runif(1),
                                        p_set_col_zero = runif(1),
                                        p_set_row_zero = runif(1)) |> add_row_names(TRUE)
    grna_matrix <- create_random_matrix(n_row = n_rows[2],
                                        n_col = n_col,
                                        p_zero = runif(1),
                                        p_set_col_zero = runif(1),
                                        p_set_row_zero = runif(1)) |> add_row_names(FALSE)
    gene_names <- generate_gene_names(gene_matrix, frac_mito = runif(1, min = 0, max = max_mito_genes))
    batch <- generate_batch(gene_matrix, sample(x = seq(1, max_n_batches), size = 1))
    curr_base_directory <- create_new_directory()
    write_sceptre_object_to_cellranger_format(gene_matrix = gene_matrix,
                                              gene_names = gene_names,
                                              grna_matrix = grna_matrix,
                                              directory = curr_base_directory,
                                              batch = batch)
    if (i %in% trials_w_vector) {
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

    return(list(base_directory = curr_base_directory, grna_matrix = grna_matrix,
                gene_matrix = gene_matrix, batch = batch, gene_names = gene_names,
                grna_target_data_frame = grna_target_data_frame))
  })

  ###############
  # run the tests
  ###############
  for (i in seq(1L, n_trials)) {
    print(paste0("Testing import from cellranger for dataset ", i))
    directories_to_load <- list.files(test_data_list[[i]]$base_directory, full.names = TRUE, pattern = "batch_*")
    # load the data from cellranger format
    if (i %in% trials_w_vector) {
      odms <- create_odm_from_cellranger(directories_to_load = directories_to_load,
                                         directory_to_write = test_data_list[[i]]$base_directory,
                                         grna_target_data_frame = test_data_list[[i]]$grna_target_data_frame,
                                         chunk_size = 10L)
    } else {
      odms <- create_odm_from_cellranger(directories_to_load = directories_to_load,
                                         directory_to_write = test_data_list[[i]]$base_directory,
                                         chunk_size = 10L)
    }
    # loop over the moalities, testing dimension and load functionality
    for (modality in c("gene", "grna")) {
      odm <- odms[[modality]]
      mem_matrix <- test_data_list[[i]][[paste0(modality, "_matrix")]]
      # if grna modality and vector supplied, perform additional processing
      vector_supplied <- i %in% trials_w_vector
      if (vector_supplied && modality == "grna") {
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
      n_row <- nrow(mem_matrix)

      # 1. check dimension
      expect_equal(dim(mem_matrix), odm@dimension)
      # 2. check feature ids
      expect_equal(rownames(mem_matrix), rownames(odm))
      # 3. check index into randomly selected rows by integer
      sample_idxs <- c(1L, sample(x = seq(1L, n_row), size = min(30, n_row), replace = FALSE), n_row)
      for (sample_idx in sample_idxs) {
        expect_equal(odm[sample_idx,], as.integer(mem_matrix[sample_idx,]))
      }
      # 4. check index into randomly selected rows by feature
      sample_features <- sample(x = rownames(mem_matrix), size = min(30, n_row), replace = FALSE)
      for (sample_feature in sample_features) {
        expect_equal(odm[sample_feature,], as.integer(mem_matrix[sample_feature,]))
      }
    }
    # check the covariates
    mem_grna_matrix <- test_data_list[[i]]$grna_matrix
    mem_gene_matrix <- test_data_list[[i]]$gene_matrix
    gene_names <- test_data_list[[i]]$gene_names
    cellwise_covariates <- odms$cellwise_covariates
    mem_batch <- test_data_list[[i]]$batch
    # i. batch
    expect_equal(as.character(cellwise_covariates$batch), as.character(mem_batch))
    # ii. n umis gene
    gene_n_umis_mem <- Matrix::colSums(mem_gene_matrix)
    gene_n_umis_disk <- cellwise_covariates$gene_n_umis
    expect_equal(as.integer(gene_n_umis_disk), as.integer(gene_n_umis_mem))
    # iii. n nonzero gene
    gene_n_nonzero_mem <- Matrix::colSums(mem_gene_matrix >= 1)
    gene_n_nonzero_disk <- cellwise_covariates$gene_n_nonzero
    expect_equal(gene_n_nonzero_mem, gene_n_nonzero_disk)
    # iv. n umis gRNA
    grna_n_umis_mem <- Matrix::colSums(mem_grna_matrix)
    grna_n_umis_disk <- cellwise_covariates$grna_n_umis
    expect_equal(as.integer(grna_n_umis_mem), as.integer(grna_n_umis_disk))
    # v. n nonzero grna
    grna_n_nonzero_mem <- Matrix::colSums(mem_grna_matrix >= 1)
    grna_n_nonzero_disk <- cellwise_covariates$grna_n_nonzero
    expect_equal(grna_n_nonzero_mem, grna_n_nonzero_disk)
    # vi. p_mito_gene
    n_umis_mito_mem <- Matrix::colSums(mem_gene_matrix[grep(pattern = "^MT-", x = gene_names),])
    p_mito_mem <- ifelse(gene_n_umis_mem == 0, 0, n_umis_mito_mem/gene_n_umis_mem)
    p_mito_disk <- cellwise_covariates$gene_p_mito
    expect_equal(p_mito_disk, p_mito_mem)
    # vii. feature_w_max_expression_grna
    grna_max_feature_res <- apply(mem_grna_matrix, 2, FUN = function(col) {
      max_feature <- which.max(col)
      max_feature_umi_count <- col[max_feature]
      list(max_feature = names(max_feature), max_feature_umi_count = max_feature_umi_count)
    })
    grna_max_feature_mem <- sapply(grna_max_feature_res, function(elem) elem$max_feature)
    grna_max_feature_disk <- cellwise_covariates$grna_feature_w_max_expression
    expect_equal(grna_max_feature_mem, grna_max_feature_disk) # ERROR
    # viii. frac max feature grna
    grna_max_exp <- sapply(grna_max_feature_res, function(elem) elem$max_feature_umi_count)
    grna_frac_max_feature_mem <- ifelse(grna_max_exp == 0, 1, grna_max_exp / grna_n_umis_mem)
    grna_frac_max_feature_disk <- cellwise_covariates$grna_frac_umis_max_feature
    expect_equal(as.numeric(grna_frac_max_feature_mem),
                 as.numeric(grna_frac_max_feature_disk))
  }
})


test_that("import data from r matrix", {
  ########################
  # define test parameters
  ########################
  n_trials <- 3L
  n_rows_range <- c(10L, 10000L)
  n_cols_range <- c(10L, 10000L)

  #################
  # create the data
  #################
  test_data_list <- lapply(seq(1L, n_trials), FUN = function(i) {
    print(paste0("Generating example dataset ", i, "."))
    n_rows <- sample(x = seq(n_rows_range[1], n_rows_range[2]), size = 2L, replace = FALSE)
    n_col <- sample(x = seq(n_cols_range[1], n_cols_range[2]), size = 1L)
    gene_matrix <- create_random_matrix(n_row = n_rows[1],
                                        n_col = n_col,
                                        p_zero = runif(1),
                                        p_set_col_zero = runif(1),
                                        p_set_row_zero = runif(1),
                                        column_access = FALSE) |> add_row_names(TRUE)
    return(gene_matrix)
  })

  ###############
  # run the tests
  ###############
  for (i in seq(1L, n_trials)) {
    print(paste0("Testing import from cellranger for dataset ", i))
    mem_matrix <- test_data_list[[i]]
    n_row <- nrow(mem_matrix)
    odm <- create_odm_from_r_matrix(mat = mem_matrix,
                                    file_to_write = paste0(tempdir(), "/gene_", i))
    # 1. check dimension
    expect_equal(dim(mem_matrix), odm@dimension)
    # 2. check feature ids
    expect_equal(rownames(mem_matrix), rownames(odm))
    # 3. check index into randomly selected rows by integer
    sample_idxs <- c(1L, sample(x = seq(1L, n_row), size = min(30, n_row), replace = FALSE), n_row)
    for (sample_idx in sample_idxs) {
      expect_equal(odm[sample_idx,], as.integer(mem_matrix[sample_idx,]))
    }
    # 4. check index into randomly selected rows by feature
    sample_features <- sample(x = rownames(mem_matrix), size = min(30, n_row), replace = FALSE)
    for (sample_feature in sample_features) {
      expect_equal(odm[sample_feature,], as.integer(mem_matrix[sample_feature,]))
    }
  }
})
