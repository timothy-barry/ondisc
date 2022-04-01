normalize_by_lib_size <- function(covariate_odm, scale_factor = 10000) {
  if (covariate_odm@post_load_function_present) stop("Data already normalized.")
  covariate_odm@post_load_function_present <- TRUE
  covariate_odm@post_load_function <- internal_normalize_by_lib_size
  return(covariate_odm)
}


internal_normalize_by_lib_size <- function(out, x, ...) {
  args <- list(...)
  # if cells were subsetted, then subset the cell libs
  cell_lib_sizes <- if ("j" %in% names(args)) x@cell_covariates$n_umis[args$j] else x@cell_covariates$n_umis
  # apply standard formula
  if (nrow(out) == 1) {
    new_out <- log(1 + out/cell_lib_sizes * scale_factor)
  } else {
    new_out <- Matrix::t(log(1 + (Matrix::t(out)/cell_lib_sizes * scale_factor)))
  }
  new_out
}


# odm_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
# metadata_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
# odm <- read_odm(odm_fp, metadata_fp)
# my_feats <- sample(get_highly_expressed_features(odm, 0.3), 100)
# covariate_odm <- odm[my_feats,]
# # add log-transformed n_umis, n_nonero
# covariate_odm <- covariate_odm %>% mutate_cell_covariates(lg_n_umis = log(n_umis), lg_n_nonzero = log(n_nonzero))
normalize_by_regression <- function(covariate_odm, covariates = c("p_mito", "lg_n_nonzero"), offset = "lg_n_umis") {
  if (covariate_odm@post_load_function_present) stop("Data already normalized.")

  # first, carry out the Poisson regressions, saving the fitted coefficients
  feature_ids <- get_feature_ids(covariate_odm)
  covariate_matrix <- get_cell_covariates(covariate_odm)
  my_form <- paste0("feature_exp ~ ", paste0(covariates, collapse = " + "), " + offset(", offset, ")")
  fitted_coefs <- lapply(X = feature_ids, function(feature_id) {
    feature_exp <- as.numeric(covariate_odm[[feature_id,]])
    fit <- stats::glm(formula = my_form, family = stats::poisson(),
                      data = dplyr::mutate(covariate_matrix, feature_exp))
    stats::coef(fit)
  })
  # save the fitted model parameters as data frame; modify column names
  fitted_coefs_df <- as.data.frame(do.call(rbind, fitted_coefs))
  fitted_coefs_df <- dplyr::rename(fitted_coefs_df, "intercept" = "(Intercept)")
  colnames(fitted_coefs_df) <- paste0("fitted_", colnames(fitted_coefs_df))

  # add the fitted_coefs_df to the feature covariate data frame; we additionally add the offset term under "misc"
  covariate_odm@feature_covariates <- dplyr::mutate(covariate_odm@feature_covariates, fitted_coefs_df)
  covariate_odm@misc <- c("offset" = offset)
}


compute_pearson_residuals <- function(out, x, ...) {
  cell_covariates <- if ("j" %in% names(args)) x@cell_covariates[args$j,] else x@cell_covariates
  # compute the fitted values

}
