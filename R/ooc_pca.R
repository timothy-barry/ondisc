#' Compute out-of-core PCA on an `ondisc` matrix
#'
#' @param covariate_odm a `covariate_ondisc_matrix` object
#' @examples
#' library(magrittr)
#' odm_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' covariate_odm <- read_odm(odm_fp, metadata_fp)
#' highly_expressed_feats <- sample(get_highly_expressed_features(odm, 0.01))
#' covariate_odm <- covariate_odm[highly_expressed_feats,]
#' highly_variable_feats <- get_highly_variable_features(covariate_odm, 2000)
#' covariate_odm <- covariate_odm[highly_variable_feats,]
#' covariate_odm <- covariate_odm %>%
#' mutate_cell_covariates(lg_n_umis = log(n_umis), lg_n_nonzero = log(n_nonzero))
#' # regress on the default covariates, i.e. p_mito and lg_n_nonzero
#' covariate_odm_norm <- normalize_by_regression(covariate_odm)
#' # run PCA covariate_odm_norm
#' compute_ooc_pca(covariate_odm_norm)
compute_ooc_pca <- function(covariate_odm) {
  MAX_N_FEATURES <- 2000
  if (nrow(covariate_odm) > MAX_N_FEATURES) {
    stop(paste0("The `covariate_odm` has ", nrow(covariate_odm)," features. However, the maximum number of features allowed in `covariate_odm` is ", MAX_N_FEATURES, "."))
  }
  # loop over cells; compute an outer product for each cell; sum these rank-1 matrices and divide by the number of cells.
  # cell_expressions %o% cell_expressions
  # perform eigen-decomposition (or spectral decomposition) on the covariance matrix (eignen function)
  # Try to determine the number of eigenvectors required to exceed a given level of variance explained (by looking at the eigenvalues)
}
