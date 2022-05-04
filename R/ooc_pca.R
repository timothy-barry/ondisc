#' Compute out-of-core PCA on an `ondisc` matrix
#'
#' @param covariate_odm_norm a `covariate_ondisc_matrix` object
#' @param fraction_variability_explained the fraction variability explained
#' @param max_num_eigenvector the maximal number of principle components
#'
#' @examples
#' library(magrittr)
#' odm_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' covariate_odm <- read_odm(odm_fp, metadata_fp)
#' highly_expressed_feats <- sample(get_highly_expressed_features(covariate_odm, 0.01))
#' covariate_odm <- covariate_odm[highly_expressed_feats,]
#' highly_variable_feats <- get_highly_variable_features(covariate_odm, 2000)
#' covariate_odm <- covariate_odm[highly_variable_feats,]
#' covariate_odm <- covariate_odm %>%
#' mutate_cell_covariates(lg_n_umis = log(n_umis), lg_n_nonzero = log(n_nonzero))
#' # regress on the default covariates, i.e. p_mito and lg_n_nonzero
#' covariate_odm_norm <- normalize_by_regression(covariate_odm)
#' # run PCA covariate_odm_norm
#' compute_ooc_pca(covariate_odm_norm, fraction_variability_explained = 0.1)
compute_ooc_pca <- function(covariate_odm_norm, fraction_variability_explained = 0.8, max_num_eigenvector = 30) {
  MAX_N_FEATURES <- 2000
  if (nrow(covariate_odm_norm) > MAX_N_FEATURES) {
    stop(paste0("The `covariate_odm_norm` has ", nrow(covariate_odm_norm)," features. However, the maximum number of features allowed in `covariate_odm` is ", MAX_N_FEATURES, "."))
  }

  # loop over cells; compute an outer product for each cell; sum these rank-1 matrices and divide by the number of cells.
  outer_products <- matrix(0, nrow = nrow(covariate_odm_norm), ncol = nrow(covariate_odm_norm))
  for (cell_idx in seq(1, ncol(covariate_odm_norm))) {
    if (cell_idx %% 100 == 0) {
      print(paste0("Processing ", cell_idx, " cells out of ", ncol(covariate_odm_norm), " cells."), quote = FALSE)
    }
    cell_expressions <- covariate_odm_norm[[,cell_idx]]
    outer_products <- outer_products + cell_expressions %o% cell_expressions
  }
  ooc_covariance <- outer_products / (ncol(covariate_odm_norm) - 1)
  eigen_decomp <- eigen(ooc_covariance)

  #Determine the number of eigenvectors required to exceed a given level of variance explained
  acc_values <- cumsum(eigen_decomp$values)
  acc_fraction <- acc_values / sum(eigen_decomp$values)
  num_vectors <- min(which(acc_fraction > fraction_variability_explained))
  if (num_vectors > max_num_eigenvector) {
    warning(paste0("Exceed maximal number of eigenvectors ", max_num_eigenvector,
                 ". Need ", num_vectors, " eigenvectors to reach fraction variability of ", round(fraction_variability_explained, 3),
                 max_num_eigenvector, " eigenvectors can only reach fraction variability of ", round(acc_fraction[max_num_eigenvector], 3)
                 ))
  }
  num_vectors <- min(num_vectors, max_num_eigenvector)

  #Compute the top k PCs: data matrix * top k eigenvectors
  top_eigenvectors <- eigen_decomp$vectors[,1:num_vectors]
  #out-of-core
  top_PCs <- matrix(0.0, nrow = ncol(covariate_odm_norm), ncol = num_vectors)
  for (cell_idx in seq(1, ncol(covariate_odm_norm))) {
    if (cell_idx %% 100 == 0) {
      print(paste0("Processing ", cell_idx, " cells out of ", ncol(covariate_odm_norm), " cells."), quote = FALSE)
    }
    cell_expressions <- t(as.numeric(covariate_odm_norm[[,cell_idx]]))
    top_PCs[cell_idx,] <- cell_expressions %*% top_eigenvectors
  }

  return(top_PCs)
}
