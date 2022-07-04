#' Get highly variable features
#'
#' Obtains the most variable features from a `covariate_ondisc_matrix`, as measured by coefficient of variation.
#'
#' @param covariate_odm a `covariate_ondisc_matrix` object
#' @param n_features the number of features to extract
#'
#' @return the top `n_features` most variable features
#' @export
#'
#' @examples
#' # Install the `ondiscdata` package before running the examples.
#' # install.packages("devtools")
#' # devtools::install_github("Katsevich-Lab/ondiscdata")
#'
#' odm_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' cov_odm <- read_odm(odm_fp, metadata_fp)
#' # first, subset to examine highly expressed features only
#' cov_odm_highly_exp <- cov_odm[get_highly_expressed_features(cov_odm),]
#' # now, extract the highly variable features
#' highly_variable_features <- get_highly_variable_features(cov_odm_highly_exp)
get_highly_variable_features <- function(covariate_odm, n_features = 250) {
  covariate_odm %>%
    get_feature_covariates() %>%
    dplyr::arrange(dplyr::desc(coef_of_variation)) %>%
    dplyr::slice(seq(1, n_features)) %>%
    row.names()
}


#' Get highly expressed features
#'
#' Obtains the set of features expressed in a least `frac_expressed` of cells.
#'
#' @param covariate_odm a `covariate_ondisc_matrix` object
#' @param frac_expressed a fraction; features expressed in this fraction (or more) of cells are returned
#'
#' @return the highly expressed features
#' @export
#'
#' @examples
#' # Install the `ondiscdata` package before running the examples.
#' # install.packages("devtools")
#' # devtools::install_github("Katsevich-Lab/ondiscdata")
#'
#' odm_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' cov_odm <- read_odm(odm_fp, metadata_fp)
#' highly_expressed_features <- get_highly_expressed_features(cov_odm)
get_highly_expressed_features <- function(covariate_odm, frac_expressed = 0.05) {
  n_cells <- ncol(covariate_odm)
  covariate_odm %>%
    get_feature_covariates() %>%
    dplyr::mutate(p_exp = n_nonzero/n_cells) %>%
    dplyr::filter(p_exp >= frac_expressed) %>%
    row.names()
}


#' Get cells with moderate sequencing depth
#'
#' Obtains the cells with a library size greater than `lower_percentile` and less than `upper_percentile`.
#'
#' @param covariate_odm a `covariate_ondisc_matrix` object
#' @param lower_percentile lower library size percentile
#' @param upper_percentile upper library size percentile
#'
#' @return a vector of cell barcodes
#' @export
#'
#' @examples
#' # Install the `ondiscdata` package before running the examples.
#' # install.packages("devtools")
#' # devtools::install_github("Katsevich-Lab/ondiscdata")
#'
#' odm_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' cov_odm <- read_odm(odm_fp, metadata_fp)
#' cells_to_keep <- get_cells_with_moderate_seq_depth(cov_odm)
get_cells_with_moderate_seq_depth <- function(covariate_odm, lower_percentile = 0.03, upper_percentile = 0.97) {
  lib_sizes <- covariate_odm %>%
    get_cell_covariates() %>%
    dplyr::pull(n_umis)
  lower_thresh <- stats::quantile(lib_sizes, lower_percentile)
  upper_thresh <- stats::quantile(lib_sizes, upper_percentile)
  covariate_odm %>%
    get_cell_covariates() %>%
    dplyr::mutate(in_range = n_umis > lower_thresh & n_umis < upper_thresh) %>%
    dplyr::filter(in_range) %>%
    row.names()
}


#' Get cells with low proportion mitochondrial reads
#'
#' Obtains the cells with proportion mitochondrial reads below a given threshold.
#'
#' @param covariate_odm a `covariate_ondisc_matrix` object
#' @param p_mito_thresh proportion mitochondrial reads threshold
#'
#' @return a vector of cell barcodes
#' @export
#'
#' @examples
#' # Install the `ondiscdata` package before running the examples.
#' # install.packages("devtools")
#' # devtools::install_github("Katsevich-Lab/ondiscdata")
#'
#' odm_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' cov_odm <- read_odm(odm_fp, metadata_fp)
#' cells_to_keep <- get_cells_with_low_p_mito(cov_odm)
get_cells_with_low_p_mito <- function(covariate_odm, p_mito_thresh = 0.1) {
  covariate_odm %>% get_cell_covariates() %>%
    dplyr::filter(p_mito <= p_mito_thresh) %>%
    row.names()
}
