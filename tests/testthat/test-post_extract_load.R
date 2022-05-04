test_that("normalize by regression", {
  for (i in seq(1, n_datasets)) {
    cov_odm <- cov_odms[[i]]
    # skip logical matrices
    if (cov_odm@ondisc_matrix@logical_mat) break()
    # extract a variable number of features
    n_feats <- sample(x = seq(1, 5L), size = 1)
    # first, get the covariate odm
    good_cells <- cov_odm %>% get_cell_covariates() %>%
      dplyr::filter(n_umis >= 1) %>%
      row.names()
    # next, get the highly expressed features
    high_exp <- get_highly_expressed_features(cov_odm, frac_expressed = 2 * 10/ncol(cov_odm))
    if (length(high_exp) < n_feats) break()
    my_features <- sample(high_exp, n_feats)
    # normalize by regression
    cov_odm <- cov_odm[my_features, good_cells]
    cov_odm <- cov_odm %>% mutate_cell_covariates(lg_n_umis = log(n_umis), lg_n_nonzero = log(n_nonzero))

    # choose covariates on which to regress
    possible_covariates <- c("p_mito", "lg_n_nonzero")
    covariates_to_regress_on <- sample(x = possible_covariates, size = sample(x = seq(0, length(possible_covariates)), size = 1))
    if (length(covariates_to_regress_on) == 0) covariates_to_regress_on <- NULL
    my_form <- paste0("feature_exp ~ ", paste0(covariates_to_regress_on, collapse = " + "), " + offset(lg_n_umis)")

    # perform the regression
    covariate_matrix <- get_cell_covariates(cov_odm)
    r_resids <- sapply(X = my_features, FUN = function(my_feat) {
      feature_exp <- as.numeric(cov_odm[[my_feat,]])
      fit <- glm(formula = my_form, family = stats::poisson(),
                 data = dplyr::mutate(covariate_matrix, feature_exp))
      pearson_resid <- resid(fit, "pearson")
    }) %>% t()
    # normalize using ondisc function
    cov_odm_normed <- normalize_by_regression(covariate_odm = cov_odm, covariates = covariates_to_regress_on)
    expect_true(all((r_resids - cov_odm_normed[[my_features,]]) < 1e-4))
    # random subset by cells
    my_cell_subset <- sample(x = seq(1, ncol(cov_odm_normed)), size = sample(x = seq(1, ncol(cov_odm_normed)), size = 1), replace = FALSE)
    expect_true(all((r_resids[,my_cell_subset] - cov_odm_normed[[my_features, my_cell_subset]]) < 1e-4))
  }
  expect_true(TRUE)
})


test_that("normalize by lib size", {
  for (i in seq(1, n_datasets)) {
    cov_odm <- cov_odms[[i]]
    scale_factor <- 10000
    good_cells <- cov_odm %>% get_cell_covariates() %>%
      dplyr::filter(n_nonzero > 0) %>% row.names()
    cov_odm <- cov_odm[,good_cells]
    cov_odm_norm <- normalize_by_lib_size(cov_odm, scale_factor)
    my_feats <- sample(x = get_feature_ids(cov_odm_norm), size = sample(x = seq(1, 10), size = 1 ))
    compare1 <- cov_odm_norm[[my_feats,]]
    compare2 <- Matrix::t(log(1 + (Matrix::t(cov_odm[[my_feats,]])/cov_odm@cell_covariates$n_umis * scale_factor)))
    expect_true(all(abs(compare1 - compare2) < 1e-4))
  }
})
