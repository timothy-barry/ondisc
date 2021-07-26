# 1. Create the matrix
n_mat <- 5
n_row_multi <- 300
n_col_multi <- sample(x = seq(100, 300), size = n_mat, replace = TRUE)
col_multi_cumsum <- c(0,cumsum(n_col_multi))
logical_mat_multi <- vector(mode = "logical", length = n_mat)
# generate the matrices using create_synthetic_data
r_mats_plus_data_multi <- vector(mode = "list", length = n_mat)
for (i in seq(1,n_mat)) {
  r_mats_plus_data_multi[[i]] <- create_synthetic_data(n_row = n_row_multi,
                                             n_col = n_col_multi[i],
                                             logical_mat = logical_mat_multi[i])
  r_mats_plus_data_multi[[i]]$features_df <- r_mats_plus_data_multi[[1]]$features_df
  r_mats_plus_data_multi[[i]]$features_fp <- r_mats_plus_data_multi[[1]]$features_fp
}
mtx_fp_multi <- sapply(X = r_mats_plus_data_multi, function(i) i$matrix_fp)
barcodes_fp_multi <- sapply(X = r_mats_plus_data_multi, function(i) i$barcodes_fp)
features_fp_multi <- r_mats_plus_data_multi[[1]]$features_fp

# 2.create ondisc matrix
cov_odm_multi <- create_ondisc_matrix_from_mtx(mtx_fp = mtx_fp_multi,
                                               barcodes_fp = barcodes_fp_multi,
                                               features_fp = features_fp_multi
                                               )

# 3.compare to ground truth
r_mats_multi <- lapply(r_mats_plus_data_multi, function(l) l$r_mat)
r_mats_parent <- do.call(cbind, r_mats_multi)
test_that("Compare ground truth matrix to on_disc_matrix from list of .mtx inputs", {
    m1 <- r_mats_parent
    on_disc_mat <- cov_odm_multi$ondisc_matrix
    m2 <- on_disc_mat[[,1:ncol(on_disc_mat)]]
    m3 <- on_disc_mat[[1:nrow(on_disc_mat),]]
    m4 <- on_disc_mat[[1:nrow(on_disc_mat),1:ncol(on_disc_mat)]]
    expect_true(all(m1 == m2))
    expect_true(all(m2 == m3))
    expect_true(all(m3 == m4))
})

# 4. check the covariates
test_that("check covariates for list of .mtx inputs", {
  r_mat <- r_mats_parent

  # first, look at cell-wise stats
  cell_covariates <-  cov_odm_multi$cell_covariates
  for (curr_col in colnames(cell_covariates)) {
    # n nonzero
    if (curr_col == "n_nonzero") {
      test <- Matrix::colSums(r_mat >= 1)
    } else if (curr_col == "n_umis") {
      test <-  Matrix::colSums(r_mat)
    } else if (curr_col == "p_mito") {
      gene_names <- get_feature_names(cov_odm_multi$ondisc_matrix)
      mt_genes <- grep(pattern = "^MT-", x = gene_names)
      mt_counts <- r_mat[mt_genes,]
      test <- Matrix::colSums(mt_counts)/Matrix::colSums(r_mat)
    }
    expect_equal(test, cell_covariates[,curr_col])
  }

  # next, look at gene-wise stats
  feature_covariates <- cov_odm_multi$feature_covariates
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
})
