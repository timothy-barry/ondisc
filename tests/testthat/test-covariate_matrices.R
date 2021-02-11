test_that("check covariates", {
  for (i in seq(1, n_datasets)) {
    r_mat <- r_mats[[i]]

    # first, look at cell-wise stats
    cell_covariates <- cov_odms[[i]]@cell_covariates
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
    feature_covariates <- cov_odms[[i]]@feature_covariates
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
