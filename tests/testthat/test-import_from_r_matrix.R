test_that("import data from R matrix", {
  ########################
  # define test parameters
  ########################
  n_trials <- 10L
  matrix_classes <- c("RsparseMatrix", "CsparseMatrix", "TsparseMatrix", "matrix")
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
                                        matrix_class = "matrix") |> add_row_names("gene")
    curr_matrix_class <- sample(matrix_classes, 1)
    gene_matrix <- as(gene_matrix, curr_matrix_class)
    return(gene_matrix)
  })

  ###############
  # run the tests
  ###############
  for (i in seq(1L, n_trials)) {
    print(paste0("Testing import from R matrix for dataset ", i))
    mem_matrix <- test_data_list[[i]]
    n_row <- nrow(mem_matrix)
    n_nonzero <- if (methods::is(mem_matrix, "matrix")) {
      sum(mem_matrix)
    } else {
      sum(mem_matrix@x)
    }

    odm <- create_odm_from_r_matrix(mat = mem_matrix,
                                    file_to_write = paste0(tempdir(), "/gene_", i, ".odm"),
                                    chunk_size = min(1000L, n_nonzero))
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
