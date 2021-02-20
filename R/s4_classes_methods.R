# Class definition and methods for ondisc_matrix


# Classes and constructors
##########################
#' `ondisc_matrix` class constructor
#'
#' An ondisc_matrix represents a feature-by-cell expression matrix stored on disk. Use this function to obtain an `ondisc_matrix` from an already-initialized on-disk ondisc_matrix.h5 file.
#'
#' @slot h5_file path to an expression matrix (created by one of the initialization functions in this package) stored on disk
#' @slot cell_subset an integer vector storing the cells currently in use
#' @slot feature_subset an integer vector storing the features currently in use
#'
#' @export ondisc_matrix
#' @examples
#' # NOTE: You must create the HDF5 file "expressions.h5" to run this example.
#' # Navigate to the help file of "create_ondisc_matrix_from_mtx"
#' # (via ?create_ondisc_matrix_from_mtx), and execute the code in the first code block.
#'
#' h5_fp <- paste0(tempdir(), "/expressions.h5")
#' if (file.exists(h5_fp)) {
#' odm <- ondisc_matrix(h5_file = h5_fp)
#' }
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


#' metadata_ondisc_matrix class constructor
#'
#' A metadata_ondisc_matrix is an object that stores and ondisc_matrix, along with its cell-specific and feature-specific covariate matrices.
#'
#' @slot ondisc_matrix an ondisc_matrix.
#' @slot cell_covariates a data frame on cell covariates.
#' @slot feature_covariates a data frame of feature covariates.
#'
#' @return a metadata_ondisc_matrix object
#' @export
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
metadata_ondisc_matrix <- setClass("metadata_ondisc_matrix",
                          slots = list(ondisc_matrix = "ondisc_matrix",
                                       cell_covariates = "data.frame",
                                       feature_covariates = "data.frame"))

#' @export
multimodal_ondisc_matrix <- setClass("multimodal_ondisc_matrix", slots = list(modalities = "list",
                                                          global_cell_covariates = "data.frame"))


#' Constructor for multimodal_ondisc_matrix
#'
#' Constructs a multimodal_ondisc_matrix from a list of metadata_ondisc_matrix objects
#'
#' @param covariate_ondisc_matrix_list a named list containing covariate_ondisc_matrices; the names are taken to be the names of the modalities
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
#'     perturbations <- readRDS(expression_fp)
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


# Basic information extraction methods
######################################

#' Return dimension
#'
#' Return dimension of ondisc_matrix.
#' @param x an ondisc_matrix
#' @export
#' @return an integer vector containing the dimension of the matrix
#' @examples
#' # NOTE: You must create the HDF5 file "expressions.h5" to run this example.
#' # Navigate to the help file of "create_ondisc_matrix_from_mtx"
#' # (via ?create_ondisc_matrix_from_mtx), and execute the code in the first code block.
#' h5_fp <- paste0(tempdir(), "/expressions.h5")
#' if (file.exists(h5_fp)) {
#' odm <- ondisc_matrix(h5_file = h5_fp)
#' dim(odm)
#' }
setMethod("dim", signature("ondisc_matrix"), function(x) get_dim(x))


#' Print basic information to the console
#'
#' @param object an object of class ondisc_matrix, covaraite_ondisc_matrix, or multimodal_ondisc_matrix
#' @return NULL
#' @export
#' @name show
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
    paste0(i, ". ", crayon::blue(modalities[i]), ": ") %>% cat()
    show(object@modalities[[i]])
  }
  cat(paste("Plus, a global cell-covariate matrix with",
             crayon::blue(ncol(object@global_cell_covariates)),
             "features."))
})


#' Print the first few rows and columns
#' @export
#' @param x an on_disc_mnatrix
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


# Subset methods and doc
########################

#' Subset an `ondisc_matrix` with `[`
#'
#' The user can pass logical, character, of numeric vectors to \code{`[`}. Character vectors correspond to feature IDs (for rows) and cell barcodes (for columns).
#' @param x An ondisc_matrix object
#' @param i Vector (numeric, logical, or character) indicating features to keep
#' @param j Vector (numeric, logical, or character) indicating cells to keep
#' @param drop not used
#' @return A subset ondisc_matrix object.
#' @name subset-odm
#' @examples
#' # NOTE: You must create the HDF5 file "expressions.h5" to run this example.
#' # Navigate to the help file of "create_ondisc_matrix_from_mtx"
#' # (via ?create_ondisc_matrix_from_mtx), and execute the code in the first code block.
#'
#' h5_fp <- paste0(tempdir(), "/expressions.h5")
#' if (file.exists(h5_fp)) {
#' odm <- ondisc_matrix(h5_file = h5_fp)
#' # keep cells 100-110
#' x <- odm[,100:110]
#' # keep the following genes:
#' x <- odm[c("ENSG00000188305", "ENSG00000257284", "ENSG00000251655"),]
#' # keep the following cells:
#' x <- odm[,c("CTTAGGACACTGGCGT-1", "AAAGGATTCACATCAG-1")]
#' # keep all genes except ENSG00000188305 and ENSG00000257284
#' x <- odm[!(get_feature_ids(odm) %in% c("ENSG00000188305", "ENSG00000257284")),]
#' }
NULL

#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "ondisc_matrix", i = "missing", j = "missing", drop = "missing"),
          definition = function(x, i, j, drop) return(x))

#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "ondisc_matrix", i = "ANY", j = "missing", drop = "missing"),
          definition = function(x, i, j) subset_by_feature_or_cell(x = x, idx = i, subset_on_cell = FALSE))

#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "ondisc_matrix", i = "missing", j = "ANY", drop = "missing"),
          definition = function(x, i, j) subset_by_feature_or_cell(x = x, idx = j, subset_on_cell = TRUE))

#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "ondisc_matrix", i = "ANY", j = "ANY", drop = "missing"),
          definition = function(x, i, j) {
            subset_by_feature_or_cell(x = x, idx = i, subset_on_cell = FALSE) %>% subset_by_feature_or_cell(x = ., idx = j, subset_on_cell = TRUE)
          })


#' Subset a `metadata_ondisc_matrix` with `[`
#'
#' The user can pass logical, character, of numeric vectors to \code{`[`}. Character vectors correspond to feature IDs (for rows) and cell barcodes (for columns).
#' @param x A metadata_ondisc_matrix object
#' @param i Vector (numeric, logical, or character) indicating features to keep
#' @param j Vector (numeric, logical, or character) indicating cells to keep
#' @param drop not used
#' @return A subset metadata_ondisc_matrix object.
#' @name subset-covariate-odm
#' @examples
#' # NOTE: You must create the HDF5 file "expressions.h5" and the RDS file
#' # "expressions.rds" to run this example. Navigate to the help file of
#' # "create_ondisc_matrix_from_mtx" (via ?create_ondisc_matrix_from_mtx),
#' # and execute the code in the first code block.
#' expressions_fp <- paste0(tempdir(), "/expressions.rds")
#' if (file.exists(expressions_fp)) {
#' expressions <- readRDS(expressions_fp)
#' # keep cells 100-110
#' x <- expressions[,100:110]
#' # keep the following genes
#' x <- expressions[c("ENSG00000188305", "ENSG00000257284", "ENSG00000251655"),]
#' }
NULL

#' @rdname subset-covariate-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "metadata_ondisc_matrix", i = "ANY", j = "ANY", drop = "missing"),
          definition = function(x, i, j, drop) {
            x@ondisc_matrix <- x@ondisc_matrix[i,j]
            x@cell_covariates <- x@cell_covariates[j,,drop=FALSE]
            x@feature_covariates <- x@feature_covariates[i,,drop=FALSE]
            return(x)
          })

#' @rdname subset-covariate-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "metadata_ondisc_matrix", i = "ANY", j = "missing", drop = "missing"),
          definition = function(x, i, j, drop) {
            x@ondisc_matrix <- x@ondisc_matrix[i,]
            x@cell_covariates <- x@feature_covariates[i,,drop=FALSE]
            return(x)
          })

#' @rdname subset-covariate-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "metadata_ondisc_matrix", i = "missing", j = "ANY", drop = "missing"),
          definition = function(x, i, j, drop) {
            x@ondisc_matrix <- x@ondisc_matrix[,j]
            x@feature_covariates <- x@cell_covariates[j,,drop=FALSE]
            return(x)
          })

#' @rdname subset-covariate-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "metadata_ondisc_matrix", i = "missing", j = "missing", drop = "missing"),
          definition = function(x, i, j, drop) return(x))


#' Subset a `multimodal_ondisc_matrix` with `[`
#'
#' The user can pass logical, character, of numeric vectors to \code{`[`}. Character vectors correspond to cell barcodes.
#' Multimodal_ondisc_matrix objects can be subset only by cell, not feature, as different modalities (in general) have different numbers of features.
#'
#' @param x A multimodal_ondisc_matrix object
#' @param i not used
#' @param j Vector (numeric, logical, or character) indicating cells to keep
#' @param drop not used
#' @return A subset metadata_ondisc_matrix object.
#' @name subset-multimodal-odm
#' @examples
#' # NOTE: You must create the RDS files "expressions.rds" and
#' # "perturbations.rds" to run this example. Navigate to the help file of
#' # "create_ondisc_matrix_from_mtx" (via ?create_ondisc_matrix_from_mtx),
#' # and execute both code blocks.
#' expression_fp <- paste0(tempdir(), "/expressions.rds")
#' perturbations_fp <- paste0(tempdir(), "/perturbations.rds")
#' if (file.exists(expression_fp) && file.exists(perturbations_fp)) {
#'     expressions <- readRDS(expression_fp)
#'     perturbations <- readRDS(expression_fp)
#'     crispr_experiment <- multimodal_ondisc_matrix(list(expressions = expressions,
#'     perturbations = perturbations))
#'     # Keep all cells except 10,100, and 105.
#'     x <- crispr_experiment[,-c(10,100,105)]
#' }
NULL

multimodal_subset_error <- function() {
  "The operation [i,] is not valid, because multimodal_ondisc_matrix objects can
  be subset only by cell, not feature. Use [,i] to subset by cell instead."
}

#' @rdname subset-multimodal-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "multimodal_ondisc_matrix", i = "missing", j = "missing", drop = "missing"),
          definition = function(x, i, j, drop) return(x))

#' @rdname subset-multimodal-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "multimodal_ondisc_matrix", i = "missing", j = "ANY", drop = "missing"),
          definition = function(x, i, j, drop) {
            n_modalities <- length(x@modalities)
            for (ctr in seq(1, n_modalities)) x@modalities[[ctr]] <- x@modalities[[ctr]][,j]
            return(x)
          })

#' @rdname subset-multimodal-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "multimodal_ondisc_matrix", i = "ANY", j = "missing", drop = "missing"),
          definition = function(x, i, j, drop) stop(multimodal_subset_error()))

#' @rdname subset-multimodal-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "multimodal_ondisc_matrix", i = "ANY", j = "ANY", drop = "missing"),
          definition = function(x, i, j, drop) stop(multimodal_subset_error()))


# Extract expression data methods
#################################

#' Pull a submatrix into memory with `[[`
#'
#' The user can pass logical, character, of numeric vectors to \code{`[[`}. Character vectors correspond to feature IDs (for rows) and cell barcodes (for columns).
#' @param x An ondisc_matrix object
#' @param i Vector (numeric, logical, or character) indicating features to keep
#' @param j Vector (numeric, logical, or character) indicating cells to keep
#' @return A matrix (of class Matrix)
#' @examples
#' # NOTE: You must create the HDF5 file "expressions.h5" to run this example.
#' # Navigate to the help file of "create_ondisc_matrix_from_mtx"
#' # (via ?create_ondisc_matrix_from_mtx), and execute the code in the first code block.
#'
#' h5_fp <- paste0(tempdir(), "/expressions.h5")
#' if (file.exists(h5_fp)) {
#' odm <- ondisc_matrix(h5_file = h5_fp)
#' # extract cells 100-110
#' x <- odm[[,100:110]]
#' # extract the following genes:
#' x <- odm[[c("ENSG00000188305", "ENSG00000257284", "ENSG00000251655"),]]
#' # extract the following cells:
#' x <- odm[[,c("CTTAGGACACTGGCGT-1", "AAAGGATTCACATCAG-1")]]
#' }
#' @name extract-odm
NULL

# 1. Extract nothing (return error).
#' @export
#' @rdname extract-odm
setMethod(f = "[[",
          signature = signature(x = "ondisc_matrix", i = "missing", j = "missing"),
          definition = function(x, i, j) stop("Specify row or column indices to extract a sub-matrix."))

# 2. Extract by feature
#' @export
#' @rdname extract-odm
setMethod(f = "[[",
          signature = signature(x = "ondisc_matrix", i = "ANY", j = "missing"),
          definition = function(x, i, j) subset_by_feature_or_cell(x = x, idx = i, subset_on_cell = FALSE) %>% extract_matrix())

# 3. Extract by cell
#' @export
#' @rdname extract-odm
setMethod(f = "[[",
          signature = signature(x = "ondisc_matrix", i = "missing", j = "ANY"),
          definition = function(x, i, j) subset_by_feature_or_cell(x = x, idx = j, subset_on_cell = TRUE) %>% extract_matrix())

# 4. Extract by both feature and cell
#' @export
#' @rdname extract-odm
setMethod(f = "[[",
          signature = signature(x = "ondisc_matrix", i = "ANY", j = "ANY"),
          definition = function(x, i, j) subset_by_feature_or_cell(x = x, idx = i, subset_on_cell = FALSE) %>% subset_by_feature_or_cell(x = ., idx = j, subset_on_cell = TRUE) %>% extract_matrix())

# covaraiate odm
################

extract_covariate_odm_error <- function() {
  "You cannot use the [[,]] operator on this object, because it is not an ondisc_matrix.
  To access the ondisc_matrix stored within this object, use the @ symbol and extract the
  field `ondisc_matrix.`"
}

# Extract covariate_odm.
#' @export
#' @rdname extract-odm
setMethod(f = "[[",
          signature = signature(x = "metadata_ondisc_matrix", i = "missing", j = "missing"),
          definition = function(x, i, j) stop(extract_covariate_odm_error()))

#' @export
#' @rdname extract-odm
setMethod(f = "[[",
          signature = signature(x = "metadata_ondisc_matrix", i = "ANY", j = "missing"),
          definition = function(x, i, j) stop(extract_covariate_odm_error()))

#' @export
#' @rdname extract-odm
setMethod(f = "[[",
          signature = signature(x = "metadata_ondisc_matrix", i = "missing", j = "ANY"),
          definition = function(x, i, j) stop(extract_covariate_odm_error()))

#' @export
#' @rdname extract-odm
setMethod(f = "[[",
          signature = signature(x = "metadata_ondisc_matrix", i = "ANY", j = "ANY"),
          definition = function(x, i, j) stop(extract_covariate_odm_error()))
