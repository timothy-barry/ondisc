# 1. mutate covariate matrix functions

#' Mutate covariate matrices.
#'
#' Functions to mutate cell-specific or feature-specific covariate matrix. Returns an `ondisc` object.
#'
#' @name mutate-covariates
#' @param x an object of class `ondisc_matrix`, `metadata_ondisc_matrix`, or `multimodal_ondisc_matrix`.
#' @param ... arguments to dplyr::mutate
#' @return an updated `ondisc` object
NULL

# Generic functions
#' @export
#' @rdname mutate-covariates
setGeneric("mutate_cell_covariates", function(x, ...) standardGeneric("mutate_cell_covariates"))

#' @export
#' @rdname mutate-covariates
setGeneric("mutate_feature_covariates", function(x, ...) standardGeneric("mutate_feature_covariates"))

# ondisc_matrix functions
#' @rdname mutate-covariates
#' @export
setMethod("mutate_cell_covariates", signature("ondisc_matrix"), function(x, ...) stop("This ondisc_matrix object does not have cell covariates."))

#' @rdname mutate-covariates
#' @export
setMethod("mutate_feature_covariates", signature("ondisc_matrix"), function(x, ...) stop("This ondisc_matrix object does not have feature covariates."))

# metadata_ondisc_matrix functions
#' @rdname mutate-covariates
#' @export
setMethod("mutate_cell_covariates", signature("metadata_ondisc_matrix"), function(x, ...) {
  x@cell_covariates <- dplyr::mutate(x@cell_covariates, ...)
  return(x)
})

#' @rdname mutate-covariates
#' @export
setMethod("mutate_feature_covariates", signature("metadata_ondisc_matrix"), function(x, ...) {
  x@feature_covariates <- dplyr::mutate(x@feature_covariates, ...)
  return(x)
})

# multimodal_ondisc_matrix functions
#' @rdname mutate-covariates
#' @export
setMethod("mutate_cell_covariates", signature("multimodal_ondisc_matrix"), function(x, ...) {
  x@global_cell_covariates <- dplyr::mutate(x@global_cell_covariates, ...)
  return(x)
})

#' @rdname mutate-covariates
#' @export
setMethod("mutate_feature_covariates", signature("multimodal_ondisc_matrix"), function(x, ...) {
  stop("mutate_feature_covariates not supported for objects of class multimodal_ondisc_matrix.")
})


# 2. Get modality
#' Get modality
#'
#' Obtains a given modality from a multimodal_ondisc_matrix object by name.
#'
#' @param multimodal_mat a multimodal_ondisc_matrix object.
#' @param modality_name name of modality to extract.
#'
#' @return the requested modality
#' @export
get_modality <- function(multimodal_mat, modality_name) {
  return(multimodal_mat@modalities[[modality_name]])
}

# 3. Get covariate matrices
#' Get covariate matrices.
#'
#' Functions to get cell-specific or feature-specific covariate matrices from `ondisc` objects.
#'
#' @name get-covariates
#' @param x an object of class `ondisc_matrix`, `metadata_ondisc_matrix`, or `multimodal_ondisc_matrix`.
#' @return a covariate matrix (in data frame form)
NULL

# Generic functions
#' @export
#' @rdname get-covariates
setGeneric("get_cell_covariates", function(x) standardGeneric("get_cell_covariates"))

#' @export
#' @rdname get-covariates
setGeneric("get_feature_covariates", function(x) standardGeneric("get_feature_covariates"))

# ondisc_matrix functions
#' @rdname get-covariates
#' @export
setMethod("get_cell_covariates", signature("ondisc_matrix"), function(x) stop("This ondisc_matrix object does not have cell covariates."))

#' @rdname get-covariates
#' @export
setMethod("get_feature_covariates", signature("ondisc_matrix"), function(x) stop("This ondisc_matrix object does not have feature covariates."))

# metadata_ondisc_matrix functions
#' @rdname get-covariates
#' @export
setMethod("get_cell_covariates", signature("metadata_ondisc_matrix"), function(x) x@cell_covariates)

#' @rdname get-covariates
#' @export
setMethod("get_feature_covariates", signature("metadata_ondisc_matrix"), function(x) x@feature_covariates)

# multimodal_ondisc_matrix functions
#' @rdname get-covariates
#' @export
setMethod("get_cell_covariates", signature("multimodal_ondisc_matrix"), function(x) x@global_cell_covariates)

#' @rdname get-covariates
#' @export
setMethod("get_feature_covariates", signature("multimodal_ondisc_matrix"), function(x) {
  stop("get_feature_covariates not supported for objects of class multimodal_ondisc_matrix.")
})
