test_that("gets, no subsets", {
  for (i in seq(1, n_datasets)) {
    on_disc_mat <- cov_odms[[i]]@ondisc_matrix
    features_df <- r_mats_plus_metadata[[i]]$features_df
    all(features_df$gene_names == get_feature_names(on_disc_mat)) %>% expect_true()
    all(features_df$gene_ids == get_feature_ids(on_disc_mat)) %>% expect_true()
    barcodes <- r_mats_plus_metadata[[i]]$barcodes
    all(barcodes == get_cell_barcodes(on_disc_mat)) %>% expect_true()
    }
})


test_that("gets after subset", {
  for (i in  seq(1, n_datasets)) {
    on_disc_mat <- cov_odms[[i]]@ondisc_matrix
      for (j in 1:n_reps) {
        subset_size_col <- sample(1:(ceiling(ncol(on_disc_mat)/10)), 1)
        subset_size_row <- sample(1:(ceiling(nrow(on_disc_mat)/10)), 1)
        col_names <- get_cell_barcodes(on_disc_mat)
        row_names <- get_feature_ids(on_disc_mat)
        # subset a first time
        t1 <- on_disc_mat[,col_names]
        expect_true(all(get_cell_barcodes(t1) == col_names))
        t2 <- on_disc_mat[row_names,]
        expect_true(all(get_feature_ids(t2) == row_names))
        # subset a second time
        subset_size_col_2 <- sample(1:length(col_names), 1)
        subset_size_row_2 <- sample(1:length(row_names), 1)
        col_names <- sample(col_names, subset_size_col_2)
        row_names <- sample(row_names, subset_size_row_2)
        t1 <- t1[,col_names]
        expect_true(all(get_cell_barcodes(t1) == col_names))
        t2 <- t2[row_names,]
        expect_true(all(get_feature_ids(t2) == row_names))
      }
    }
})


test_that("covariate matrix gets", {
  on_disc_mat <- cov_odms[[1]]@ondisc_matrix
  expect_error(get_feature_covariates(on_disc_mat))
  expect_error(get_cell_covariates(on_disc_mat))
  expect_error(mutate_cell_covariates(on_disc_mat))
  expect_error(mutate_feature_covariates(on_disc_mat))
})
