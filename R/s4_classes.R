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
                                        underlying_dimension = "integer",
                                        cell_subset = "integer",
                                        feature_subset = "integer",
                                        feature_ids = "character",
                                        feature_names = "character",
                                        cell_barcodes = "character",
                                        odm_id = "integer"),
                           prototype = list(h5_file = NA_character_,
                                            logical_mat = FALSE,
                                            underlying_dimension = NA_integer_,
                                            cell_subset = NA_integer_,
                                            feature_subset = NA_integer_,
                                            feature_ids = NA_character_,
                                            feature_names = NA_character_,
                                            cell_barcodes = NA_character_,
                                            odm_id = NA_integer_)
                          )

###########################
# 2. covariate_ondisc_matrix
###########################

#' `covariate_ondisc_matrix` class
#'
#' A `covariate_ondisc_matrix` stores an `ondisc_matrix`, along with cell-specific and feature-specific covariate matrices.
#'
#' @slot ondisc_matrix an ondisc_matrix.
#' @slot cell_covariates a data frame of cell covariates.
#' @slot feature_covariates a data frame of feature covariates.
covariate_ondisc_matrix <- setClass("covariate_ondisc_matrix",
                          slots = list(ondisc_matrix = "ondisc_matrix",
                                       cell_covariates = "data.frame",
                                       feature_covariates = "data.frame"))


#' `covariate_ondisc_matrix` constructor
#'
#' Construct a `covariate_ondisc_matrix` by passing an `ondisc_matrix`, along with its associated `cell_covariates` and `feature_covariates`.
#'
#' @param ondisc_matrix an `ondisc_matrix`.
#' @param cell_covariates a data frame storing the cell-specific covariates.
#' @param feature_covariates a data frame storing the feature-specific covariates.
#'
#' @return a `covariate_ondisc_matrix`.
#' @export
covariate_ondisc_matrix <- function(ondisc_matrix, cell_covariates, feature_covariates) {
  out <- new("covariate_ondisc_matrix")
  out@ondisc_matrix <- ondisc_matrix
  out@cell_covariates <- cell_covariates
  out@feature_covariates <- feature_covariates
  return(out)
}


#' `multimodal_ondisc_matrix` class
#'
#' A `multimodal_ondisc_matrix` represents multimodal data.
#'
#' @slot modalities a list containing `covariate_ondisc_matrix` objects representing different modalities.
#' @slot global_cell_covariates a data frame containing the cell-specific covariates pooled across all modalities.
multimodal_ondisc_matrix <- setClass("multimodal_ondisc_matrix", slots = list(modalities = "list",
                                                          global_cell_covariates = "data.frame"))

#' `multimodal_ondisc_matrix` constructor
#'
#' Construct a `multimodal_ondisc_matrix` from a list of `covariate_ondisc_matrix` objects.
#'
#' @param covariate_ondisc_matrix_list a named list containing `metadata_ondisc_matrices`; the names are taken to be the names of the modalities.
#'
#' @return a multimodal_ondisc_matrix
#' @export
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
