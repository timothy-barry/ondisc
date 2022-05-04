test_that("check ooc pca", {
  odm_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
  metadata_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
  covariate_odm <- read_odm(odm_fp, metadata_fp)
  highly_expressed_feats <- sample(get_highly_expressed_features(covariate_odm, 0.01))
  covariate_odm <- covariate_odm[highly_expressed_feats,]
  highly_variable_feats <- get_highly_variable_features(covariate_odm, 500)
  covariate_odm <- covariate_odm[highly_variable_feats,]
  covariate_odm <- covariate_odm %>%
    mutate_cell_covariates(lg_n_umis = log(n_umis), lg_n_nonzero = log(n_nonzero))
  # regress on the default covariates, i.e. p_mito and lg_n_nonzero
  covariate_odm_norm <- normalize_by_regression(covariate_odm)

  # first, compute ooc pca
  ooc_PCs <- compute_ooc_pca(covariate_odm_norm, fraction_variability_explained = 0.1)

  # next, run in-core pca
  data_matrix <- covariate_odm_norm[[1:nrow(covariate_odm_norm), 1:ncol(covariate_odm_norm)]]
  data_matrix <- as.matrix(t(data_matrix))
  in_core_covariance <- cov(data_matrix)
  eigen_vectors <- eigen(in_core_covariance)
  top_eigenvectors <- eigen_vectors$vectors[, 1:ncol(ooc_PCs)]

  in_core_PCs <- data_matrix %*% top_eigenvectors

  ooc_svd <- svd(ooc_PCs)
  in_core_svd <- svd(in_core_PCs)
  expect_true(abs(ooc_svd$d[1] - in_core_svd$d[1]) < 1e-5)
})
