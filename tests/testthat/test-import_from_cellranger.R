test_that("import data from cellranger", {
  ########################
  # define test parameters
  ########################
  n_trials <- 3L
  n_rows_range <- c(10L, 10000L)
  n_cols_range <- c(10L, 10000L)
  max_mito_genes <- 0.1
  max_n_batches <- 5

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
    return(list(base_directory = curr_base_directory, grna_matrix = grna_matrix,
                gene_matrix = gene_matrix, batch = batch, gene_names = gene_names))
  })

  ###############
  # run the tests
  ###############
  for (i in seq(1L, n_trials)) {
    print(paste0("Testing import from cellranger for dataset ", i))
    directories_to_load <- list.files(test_data_list[[i]]$base_directory, full.names = TRUE, pattern = "batch_*")
    # load the data from cellranger format
    odms <- create_odm_from_cellranger(directories_to_load = directories_to_load,
                                       directory_to_write = test_data_list[[i]]$base_directory)
    # loop over the moalities, testing dimension and load functionality
    for (modality in c("gene", "grna")) {
      odm <- odms[[modality]]
      mem_matrix <- test_data_list[[i]][[paste0(modality, "_matrix")]]
      n_row <- nrow(mem_matrix)
      # 1. check dimension
      expect_equal(dim(mem_matrix), odm@dimension)
      # 2. check feature ids
      expect_equal(rownames(mem_matrix), get_feature_ids(odm))
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
    n_umis_gene_mem <- Matrix::colSums(mem_gene_matrix)
    n_umis_gene_disk <- cellwise_covariates$n_umis_gene
    expect_equal(as.integer(n_umis_gene_disk), as.integer(n_umis_gene_mem))
    # iii. n nonzero gene
    n_nonzero_gene_mem <- Matrix::colSums(mem_gene_matrix >= 1)
    n_nonzero_gene_disk <- cellwise_covariates$n_nonzero_features_gene
    expect_equal(n_nonzero_gene_mem, n_nonzero_gene_disk)
    # iv. n umis gRNA
    n_umis_grna_mem <- Matrix::colSums(mem_grna_matrix)
    n_umis_grna_disk <- cellwise_covariates$n_umis_grna
    expect_equal(as.integer(n_umis_grna_mem), as.integer(n_umis_grna_disk))
    # v. n nonzero grna
    n_nonzero_grna_mem <- Matrix::colSums(mem_grna_matrix >= 1)
    n_nonzero_grna_disk <- cellwise_covariates$n_nonzero_features_grna
    expect_equal(n_nonzero_grna_mem, n_nonzero_grna_disk)
    # vi. p_mito_gene
    n_umis_mito_mem <- Matrix::colSums(mem_gene_matrix[grep(pattern = "^MT-", x = gene_names),])
    p_mito_mem <- ifelse(n_umis_gene_mem == 0, 0, n_umis_mito_mem/n_umis_gene_mem)
    p_mito_disk <- cellwise_covariates$p_mito_gene
    expect_equal(p_mito_disk, p_mito_mem)
    # vii. feature_w_max_expression_grna
    grna_max_feature_res <- apply(mem_grna_matrix, 2, FUN = function(col) {
      max_feature <- which.max(col)
      max_feature_umi_count <- col[max_feature]
      c(max_feature = max_feature - 1L, max_feature_umi_count = max_feature_umi_count)
    })
    grna_max_feature_mem <- grna_max_feature_res[1,]
    grna_max_feature_disk <- cellwise_covariates$feature_w_max_expression_grna
    expect_equal(grna_max_feature_mem, grna_max_feature_disk)
    # viii. frac max feature grna
    grna_frac_max_feature_mem <- ifelse(n_umis_grna_mem == 0, 1, grna_max_feature_res[2,] / n_umis_grna_mem)
    grna_frac_max_feature_disk <- cellwise_covariates$frac_umis_max_feature_grna
    expect_equal(grna_frac_max_feature_mem, grna_frac_max_feature_disk)
  }
})
