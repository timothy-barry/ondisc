# Class definition and methods for ondisc_matrix

##########################
# Classes and constructors
##########################

##################
# 1. ondisc_matrix
##################

#' `ondisc_matrix` class
#'
#' An `ondisc_matrix` represents a feature-by-cell expression matrix stored on-disk.
#'
#' It is best to avoid interacting with the slots of an `ondisc_matrix` directly. Instead, use the functions
#' and operators provided by the package.
#'
#' @slot h5_file path to an initialized .h5 file stored on-disk.
#' @slot logical_mat logical value indicating whether the matrix is logical.
#' @slot cell_subset integer vector recording the cells currently in use.
#' @slot feature_subset integer vector recording the features currently in use.
#' @slot underlying_dimension the dimension of the (unsubset) expression matrix.
ondisc_matrix <- setClass("ondisc_matrix",
                           slots = list(h5_file = "character",
                                        logical_mat = "logical",
                                        cell_subset = "integer",
                                        feature_subset = "integer",
                                        underlying_dimension = "integer"),
                           prototype = list(h5_file = NA_character_,
                                            logical_mat = FALSE,
                                            cell_subset = NA_integer_,
                                            feature_subset = NA_integer_,
                                            underlying_dimension = NA_integer_))

#' `ondisc_matrix` constructor
#'
#' Construct an `ondisc_matrix` from an initialized .h5 file.
#'
#' @param h5_file a .h5 file storing the on-disk portion of an initialized `ondisc_matrix` object.
#'
#' @return an initialized `ondisc_matrix` object.
#' @export
#'
#' @examples
#' # NOTE: You must create the HDF5 file "expressions.h5" to run this example.
#' # Navigate to the help file of "create_ondisc_matrix_from_mtx"
#' # (via ?create_ondisc_matrix_from_mtx), and execute the code in the first code block.
#' h5_fp <- paste0(tempdir(), "/expressions.h5")
#' if (file.exists(h5_fp)) {
#'   odm <- ondisc_matrix(h5_file = h5_fp)
#' }
ondisc_matrix <- function(h5_file) {
  out <- new(Class = "ondisc_matrix")
  out@h5_file <- h5_file
  out@underlying_dimension <- as.integer(rhdf5::h5read(file = h5_file, name = "dimension"))
  out@logical_mat <- as.logical(rhdf5::h5read(file = h5_file, name = "logical_mat"))
  return(out)
}

###########################
# 2. metadata_ondisc_matrix
###########################

#' `metadata_ondisc_matrix` class
#'
#' A `metadata_ondisc_matrix` stores an `ondisc_matrix`, along with cell-specific and feature-specific covariate matrices.
#'
#' @slot ondisc_matrix an ondisc_matrix.
#' @slot cell_covariates a data frame of cell covariates.
#' @slot feature_covariates a data frame of feature covariates.
metadata_ondisc_matrix <- setClass("metadata_ondisc_matrix",
                          slots = list(ondisc_matrix = "ondisc_matrix",
                                       cell_covariates = "data.frame",
                                       feature_covariates = "data.frame"))


#' `metadata_ondisc_matrix` constructor
#'
#' Construct a `metadata_ondisc_matrix` by passing an `ondisc_matrix`, along with its associated `cell_covariates` and `feature_covariates`.
#'
#' @param ondisc_matrix an `ondisc_matrix`.
#' @param cell_covariates a data frame storing the cell-specific covariates.
#' @param feature_covariates a data frame storing the feature-specific covariates.
#'
#' @return a `metadata_ondisc_matrix`.
#' @export
#'
#' @examples
#' # NOTE: You must create the HDF5 file "expressions.h5" and the RDS file
#' # "expressions.rds" to run this example. Navigate to the help file of
#' # "create_ondisc_matrix_from_mtx" (via ?create_ondisc_matrix_from_mtx),
#' # and execute the code in the first code block.
#' covariates_fp <- paste0(tempdir(), "/expressions.rds")
#' h5_fp <-  paste0(tempdir(), "/expressions.h5")
#' if (file.exists(covariates_fp) && file.exists(h5_fp)) {
#' covariate_odm <- readRDS(covariates_fp)
#' cell_covariate_matrix <- covariate_odm@cell_covariates
#' feature_covariate_matrix <- covariate_odm@feature_covariates
#' covariate_odm_copy <- metadata_ondisc_matrix(
#' ondisc_matrix = ondisc_matrix(h5_file = h5_fp),
#' cell_covariates = cell_covariate_matrix,
#' feature_covariates = feature_covariate_matrix)
#' }
metadata_ondisc_matrix <- function(ondisc_matrix, cell_covariates, feature_covariates) {
  out <- new("metadata_ondisc_matrix")
  out@ondisc_matrix <- ondisc_matrix
  out@cell_covariates <- cell_covariates
  out@feature_covariates <- feature_covariates
  return(out)
}


#' `multimodal_ondisc_matrix` class
#'
#' A `multimodal_ondisc_matrix` represents multimodal data.
#'
#' @slot modalities a list containing `metadata_ondisc_matrix` objects representing different modalities.
#' @slot global_cell_covariates a data frame containing the cell-specific covariates pooled across all modalities.
multimodal_ondisc_matrix <- setClass("multimodal_ondisc_matrix", slots = list(modalities = "list",
                                                          global_cell_covariates = "data.frame"))

#' `multimodal_ondisc_matrix` constructor
#'
#' Construct a `multimodal_ondisc_matrix` from a list of `metadata_ondisc_matrix` objects.
#'
#' @param covariate_ondisc_matrix_list a named list containing `metadata_ondisc_matrices`; the names are taken to be the names of the modalities.
#'
#' @return a multimodal_ondisc_matrix
#' @export
#' @examples
#' # NOTE: You must create the RDS files "expressions.rds" and
#' # "perturbations.rds" to run this example. Navigate to the help file of
#' # "create_ondisc_matrix_from_mtx" (via ?create_ondisc_matrix_from_mtx),
#' # and execute both code blocks.
#' expression_fp <- paste0(tempdir(), "/expressions.rds")
#' perturbations_fp <- paste0(tempdir(), "/perturbations.rds")
#' if (file.exists(expression_fp) && file.exists(perturbations_fp)) {
#'     expressions <- readRDS(expression_fp)
#'     perturbations <- readRDS(perturbations_fp)
#'     crispr_experiment <- multimodal_ondisc_matrix(list(expressions = expressions,
#'     perturbations = perturbations))
#' }
multimodal_ondisc_matrix <- function(covariate_ondisc_matrix_list) {
  out <- new(Class = "multimodal_ondisc_matrix")
  out@modalities <- covariate_ondisc_matrix_list
  df_list <- lapply(X = covariate_ondisc_matrix_list,
                    FUN = function(cov_odm) cov_odm@cell_covariates)
  modality_names <- names(covariate_ondisc_matrix_list)
  global_df <- combine_multimodal_dataframes(df_list, modality_names)
  out@global_cell_covariates <- global_df
  return(out)
}
