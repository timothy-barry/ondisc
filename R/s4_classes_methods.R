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

######################################
# Basic information extraction methods
######################################

#' Return dimension
#'
#' Return the dimension of an `ondisc_matrix`, `metadata_ondisc_matrix`, or `multimodal_ondisc_matrix`.
#' @param x an `ondisc_matrix`, `metadata_ondisc_matrix`, or `multimodal_ondisc_matrix`.
#' @name dim
#' @return If `x` is an `ondisc_matrix` or `metadata_ondisc_matrix`, length-two integer vector containing
#' the dimension of `x`; if `x` is a `multimodal_ondisc_matrix`, a list of integer vectors containing the dimensions
#' of the constituent modalities of `x`.
#' @examples
#' # NOTE: You must create the RDS files "expressions.rds" and
#' # "perturbations.rds" to run this example. Navigate to the help file of
#' # "create_ondisc_matrix_from_mtx" (via ?create_ondisc_matrix_from_mtx),
#' # and execute both code blocks.
#'
#' # dimension of an ondisc_matrix
#' h5_fp <- paste0(tempdir(), "/expressions.h5")
#' if (file.exists(h5_fp)) {
#' odm <- ondisc_matrix(h5_file = h5_fp)
#' dim(odm)
#' }
#'
#' # dimension of a metadata_ondic_matrix
#' expressions_fp <- paste0(tempdir(), "/expressions.rds")
#' if (file.exists(expressions_fp)) {
#' expressions <- readRDS(expressions_fp)
#' dim(expressions)
#' }
#'
#' # dimension of a multimodal_ondisc_matrix
#' expression_fp <- paste0(tempdir(), "/expressions.rds")
#' perturbations_fp <- paste0(tempdir(), "/perturbations.rds")
#' if (file.exists(expression_fp) && file.exists(perturbations_fp)) {
#'     crispr_experiment <- multimodal_ondisc_matrix(list(expressions = readRDS(expression_fp),
#'     perturbations = readRDS(perturbations_fp)))
#'     dim(crispr_experiment)
#' }
NULL

#' @export
#' @rdname dim
setMethod("dim", signature("ondisc_matrix"), function(x) get_dim(x))

#' @export
#' @rdname dim
setMethod("dim", signature("metadata_ondisc_matrix"), function(x) get_dim(x@ondisc_matrix))

#' @export
#' @rdname dim
setMethod("dim", signature("multimodal_ondisc_matrix"), function(x) lapply(x@modalities, dim))


#' Print basic information to the console
#'
#' @param object an object of class ondisc_matrix, covaraite_ondisc_matrix, or multimodal_ondisc_matrix
#' @name show
#' @return NULL
#' @examples
#' # NOTE: You must create the HDF5 file "expressions.h5" to run this example.
#' # Navigate to the help file of "create_ondisc_matrix_from_mtx"
#' # (via ?create_ondisc_matrix_from_mtx), and execute the code in the first code block.
#'
#' h5_fp <- paste0(tempdir(), "/expressions.h5")
#' if (file.exists(h5_fp)) {
#' odm <- ondisc_matrix(h5_file = h5_fp)
#' show(odm)
#' }
NULL

#' @rdname show
#' @export
setMethod("show", signature = signature("ondisc_matrix"), function(object) {
  x_dim <- dim(object)
  cat(paste0("An ondisc_matrix with ", crayon::blue(x_dim[1]), " feature", if (x_dim[1] == 1) "" else "s" ," and ", crayon::blue(x_dim[2]), " cell", if (x_dim[2] == 1) "" else "s", ".\n"))
})

#' @rdname show
#' @export
setMethod("show", signature = signature("metadata_ondisc_matrix"), function(object) {
  cell_covariates <- colnames(object@cell_covariates)
  feature_covariates <- colnames(object@feature_covariates)
  cat("A metadata_ondisc_matrix with the following components:\n")
  cat("\t"); show(object@ondisc_matrix)
  paste0("\tA cell covariate matrix with columns ", paste(crayon::blue(cell_covariates), collapse = ", "), ".\n") %>% cat()
  paste0("\tA feature covariate matrix with columns ", paste(crayon::blue(feature_covariates), collapse = ", "), ".\n") %>% cat()
})

#' @rdname show
#' @export
setMethod("show", signature = signature("multimodal_ondisc_matrix"), function(object) {
  modalities <- names(object@modalities)
  cat("A multimodal_ondisc_matrix with the following modalities:\n")
  for (i in seq(1, length(modalities))) {
    paste0(i, ". ", crayon::blue(modalities[i]), " ") %>% cat()
    show(object@modalities[[i]])
  }
  cat(paste("Plus, a global cell-covariate matrix containing",
            crayon::blue(ncol(object@global_cell_covariates)), "covariates and",
            crayon::blue(nrow(object@global_cell_covariates)), "cells."))
})


#' Print the first few rows and columns
#' @export
#' @param x an ondisc_matrix
#' @return NULL
#' @examples
#' # NOTE: You must create the HDF5 file "expressions.h5" to run this example.
#' # Navigate to the help file of "create_ondisc_matrix_from_mtx"
#' # (via ?create_ondisc_matrix_from_mtx), and execute the code in the first code block.
#'
#' h5_fp <- paste0(tempdir(), "/expressions.h5")
#' if (file.exists(h5_fp)) {
#' odm <- ondisc_matrix(h5_file = h5_fp)
#' head(odm)
#' }
setMethod("head", signature = signature("ondisc_matrix"), function(x) {
  x_dim <- dim(x)
  n_row_to_show <- min(5, x_dim[1])
  n_col_to_show <- min(6, x_dim[2])
  cat(paste0("Showing ", crayon::blue(n_row_to_show), " of ", crayon::blue(x_dim[1]), " features", if (x_dim[1] == 1) "" else "s" ," and ", crayon::blue(n_col_to_show), " of ", crayon::blue(x_dim[2])," cell", if (n_col_to_show == 1) "" else "s", ":\n"))
  print(as.matrix(x[[1:n_row_to_show, 1:n_col_to_show]]))
})
