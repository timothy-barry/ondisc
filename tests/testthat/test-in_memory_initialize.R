# Initialize metadata odms from in-memory R objects in different directories
cov_odms_from_memory <- lapply(r_mats_plus_metadata, function(l) {
  file_dir <- create_new_directory()
  # randomly choose the class of matrix
if (is.logical(r_matrix@x)) {
    r_matrix <- as(r_matrix, "lgTMatrix")
  } else {
    r_matrix <- as.matrix(l$r_mat)
    matrix_class <- sample(c("dgTMatrix", "dgRMatrix", "dgCMatrix", "matrix"), 1)
    r_matrix <- as(r_matrix, matrix_class)
  }

  metadata_odm <- create_ondisc_matrix_from_R_matrix(r_matrix = r_matrix,
                                                     barcodes = l$barcodes,
                                                     features_df = l$features_df,
                                                     odm_fp = file_dir)
})


test_that("Compare ground truth matrix to on_disc_matrix initialized in memory", {
  for (i in seq(1, n_datasets)) {
    m1 <- r_mats[[i]]
    on_disc_mat <- cov_odms_from_memory[[i]]@ondisc_matrix
    m2 <- on_disc_mat[[,1:ncol(on_disc_mat)]]
    m3 <- on_disc_mat[[1:nrow(on_disc_mat),]]
    m4 <- on_disc_mat[[1:nrow(on_disc_mat),1:ncol(on_disc_mat)]]
    expect_true(all(m1 == m2))
    expect_true(all(m2 == m3))
    expect_true(all(m3 == m4))
  }
})


test_that("Check dimension", {
  for (i in seq(1, n_datasets)) {
    original_dim <- dim(r_mats[[i]])
    test_dim <- dim(cov_odms_from_memory[[i]]@ondisc_matrix)
    expect_true(all(original_dim == test_dim))
  }
})


test_that("check covariates", {
for (i in seq(1, n_datasets)) {
  r_mat <- r_mats[[i]]

  # first, look at cell-wise stats
  cell_covariates <- cov_odms_from_memory[[i]]@cell_covariates
  for (curr_col in colnames(cell_covariates)) {
    # n nonzero
    if (curr_col == "n_nonzero") {
      test <- Matrix::colSums(r_mat >= 1)
    } else if (curr_col == "n_umis") {
      test <-  Matrix::colSums(r_mat)
    } else if (curr_col == "p_mito") {
      gene_names <- get_feature_names(cov_odms[[i]]@ondisc_matrix)
      mt_genes <- grep(pattern = "^MT-", x = gene_names)
      mt_counts <- r_mat[mt_genes,]
      test <- Matrix::colSums(mt_counts)/Matrix::colSums(r_mat)
    }
    expect_equal(test, cell_covariates[,curr_col])
  }

  # next, look at gene-wise stats
  feature_covariates <- cov_odms_from_memory[[i]]@feature_covariates
  for (curr_col in colnames(feature_covariates)) {
    if (curr_col == "n_nonzero") {
      test <- Matrix::rowSums(r_mat >= 1)
    } else if (curr_col == "mean_expression") {
      test <- Matrix::rowMeans(r_mat)
    } else if (curr_col == "coef_of_variation") {
      n_cells <- ncol(r_mat)
      my_vars <- Matrix::rowSums(r_mat^2)/n_cells - (Matrix::rowSums(r_mat)/n_cells)^2
      my_means <- Matrix::rowMeans(r_mat)
      test <- sqrt(my_vars)/my_means
    }
    expect_equal(test, feature_covariates[,curr_col])
  }
}
})

