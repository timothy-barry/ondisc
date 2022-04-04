normalize_by_lib_size <- function(covariate_odm, scale_factor = 10000) {
  if (covariate_odm@post_load_function_present) stop("Data already normalized.")
  covariate_odm@post_load_function_present <- TRUE
  covariate_odm@post_load_function <- internal_normalize_by_lib_size
  covariate_odm@misc <- list(scale_factor = scale_factor)
  return(covariate_odm)
}


internal_normalize_by_lib_size <- function(out, x, ...) {
  scale_factor <- x@misc$scale_factor
  args <- list(...)
  # if cells were subsetted, then subset the cell libs
  cell_lib_sizes <- if ("j" %in% names(args)) x@cell_covariates$n_umis[args$j] else x@cell_covariates$n_umis
  # apply standard formula
  #if (nrow(out) == 1) {
  #  new_out <- log(1 + out/cell_lib_sizes * scale_factor)
  #} else {
  #}
  Matrix::t(log(1 + (Matrix::t(out)/cell_lib_sizes * scale_factor)))
}


# odm_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
# metadata_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
# odm <- read_odm(odm_fp, metadata_fp)
# my_feats <- sample(get_highly_expressed_features(odm, 0.3), 100)
# covariate_odm <- odm[my_feats,]
# # add log-transformed n_umis, n_nonero
# covariate_odm <- covariate_odm %>% mutate_cell_covariates(lg_n_umis = log(n_umis), lg_n_nonzero = log(n_nonzero))
# covariate_odm_norm <- normalize_by_regression(covariate_odm)
normalize_by_regression <- function(covariate_odm, covariates = c("p_mito", "lg_n_nonzero"), offset = "lg_n_umis") {
  EXPRESSION_CUTOFF <- 10
  if (covariate_odm@post_load_function_present) stop("Data already normalized.")
  covariate_matrix <- get_cell_covariates(covariate_odm)

  # verify that the covariates and offset are in fact columns of the cell covariate matrix
  if (!all(c(covariates, offset) %in% colnames(covariate_matrix))) {
    stop("Ensure that `covariates` and `offset` are columns of the cell covariate matrix of `covariate_odm` (obtainable via `get_cell_covariates(covariate_odm)` and modifiable via `mutate_cell_covariates(covariate_odm)`).")
  }

  # Carry out the Poisson regressions, saving the fitted coefficients
  feature_ids <- get_feature_ids(covariate_odm)
  my_form <- paste0("feature_exp ~ ", paste0(covariates, collapse = " + "), " + offset(", offset, ")")
  fitted_coefs <- lapply(X = feature_ids, function(feature_id) {
    feature_exp <- as.numeric(covariate_odm[[feature_id,]])
    if (sum(feature_exp) <= EXPRESSION_CUTOFF) stop(paste0("The feature ", feature_id, " is too lowly expressed to normalize via regression; remove this feature from the `covariate_odm` and try again."))
    fit <- stats::glm(formula = my_form, family = stats::poisson(),
                      data = dplyr::mutate(covariate_matrix, feature_exp))
    stats::coef(fit)
  })
  # save the fitted model parameters as data frame; modify column names
  fitted_coefs_df <- as.data.frame(do.call(rbind, fitted_coefs)) %>%
    dplyr::rename(".intercept" = "(Intercept)") %>%
    dplyr::rename_with(function(col_name) paste0(".fit_", col_name), -.intercept)

  # add the fitted_coefs_df to the feature covariate data frame; we additionally add the offset term under "misc"
  covariate_odm@feature_covariates <- dplyr::mutate(covariate_odm@feature_covariates, fitted_coefs_df)
  covariate_odm@misc <- list(offset = offset, covariates = covariates)

  # finally, update fields of
  covariate_odm@post_load_function_present <- TRUE
  covariate_odm@post_load_function <- compute_pearson_residuals
  return(covariate_odm)
}


compute_pearson_residuals <- function(out, x, ...) {
  args <- list(...)
  cell_covariates <- if ("j" %in% names(args)) x@cell_covariates[args$j,] else x@cell_covariates
  feature_covariates <- if ("i" %in% names(args)) x@feature_covariates[args$i,] else x@feature_covariates
  # compute the fitted values (apply over the rows of feature covariates)
  offset <- x@misc$offset
  covariates <- x@misc$covariates
  covariates_matrix <- as.matrix(cell_covariates[,covariates])
  coefs_matrix <- as.matrix(feature_covariates[,paste0(".fit_", covariates)])
  intercept_numeric <- feature_covariates[, ".intercept",]
  offset_numeric <- cell_covariates[[offset]]

  # compute Pearson residuals
  #if (nrow(out) == 1) {
  #  mu_hat <- exp(intercept_numeric + as.numeric(covariates_matrix %*%  t(coefs_matrix)) + offset_numeric)
  #  ret <- (out - mu_hat)/sqrt(mu_hat)
  #} else {
  #}

  mu_hat <- exp(intercept_numeric + t(covariates_matrix %*% t(coefs_matrix) + offset_numeric))
  ret <- (out - mu_hat)/sqrt(mu_hat)
  return(ret)
}
