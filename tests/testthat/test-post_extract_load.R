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
    cov_odm_normed <- normalize_by_regression(covariate_odm = cov_odm)

    # run regression on the selected features and manually extract Pearson resids
    my_form <- paste0("feature_exp ~ ", paste0(cov_odm_normed@misc$covariates, collapse = " + "), " + offset(", cov_odm_normed@misc$offset, ")")
    covariate_matrix <- get_cell_covariates(cov_odm)

    r_resids <- sapply(X = my_features, FUN = function(my_feat) {
      feature_exp <- as.numeric(cov_odm[[my_feat,]])
      fit <- glm(formula = my_form, family = stats::poisson(),
                 data = dplyr::mutate(covariate_matrix, feature_exp))
      pearson_resid <- resid(fit, "pearson")
    }) %>% t()
    expect_true(all((r_resids - cov_odm_normed[[my_features,]]) < 1e-4))
    # random subset by cells
    my_cell_subset <- sample(x = seq(1, ncol(cov_odm_normed)), size = sample(x = seq(1, ncol(cov_odm_normed)), size = 1), replace = FALSE)
    expect_true(all((r_resids[,my_cell_subset] - cov_odm_normed[[my_features, my_cell_subset]]) < 1e-4))
  }
  expect_true(TRUE)
})
