# 1. mutate covariate matrix functions

#' Mutate covariate matrices.
#'
#' Functions to mutate cell-specific or feature-specific covariate matrices in `covariate_ondisc_matrix` or `multimodal_ondisc_matrix` objects.
#' For `covariate_ondisc_matrix`, `mutate_feature_covariates` can be used to modify the feature-specific covariate matrix and `mutate_cell_covariates` can be used to modify cell-specific covariate matrix.
#' An updated `covariate_ondisc_matrix` will be returned.
#'
#' For `multimodal_ondisc_matrix`, `mutate_cell_covariates` can be used to modify cell-specific covariate matrix but `mutate_feature_covariates` not supported is not supported.
#' The modification is made on the `global_cell_covariates` but not the cell-specific or feature-specific covariate matrices in each respective modal.
#' An updated `multimodal_ondisc_matrix` will be returned.
#'
#' @name mutate-covariates
#' @param x an object of class `covariate_ondisc_matrix` or `multimodal_ondisc_matrix`.
#' @param ... arguments to dplyr::mutate
#' @return an updated `covariate_ondisc_matrix` or `multimodal_ondisc_matrix` object as the same type as the input
#' @examples
#' # Install the `ondiscdata` package to run the examples.
#' # devtools::install_github("timothy-barry/ondiscdata")
#'
#' ####################################
#' # EXAMPLE 1: covariate_ondisc_matrix
#' ####################################
#' odm_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' odm <- read_odm(odm_fp, metadata_fp)
#'
#' # add one column for cell covariate
#' odm_p_feature_expressed <- mutate_cell_covariates(odm, p_feature_expressed = n_nonzero/nrow(odm))
#' # delete coef_of_variation from feature covariate
#' odm_no_coef_of_variation <- mutate_feature_covariates(odm, coef_of_variation = NULL)
#'
#' #####################################
#' # EXAMPLE 2: multimodal_ondisc_matrix
#' #####################################
#' odm_gene_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_gene_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' odm_gene <- read_odm(odm_gene_fp, metadata_gene_fp)
#' odm_grna_fp <- system.file("extdata",
#' "odm/grna_assignment/matrix.odm", package = "ondiscdata")
#' metadata_grna_fp <- system.file("extdata",
#' "odm/grna_assignment/metadata.rds", package = "ondiscdata")
#' odm_grna <- read_odm(odm_grna_fp, metadata_grna_fp)
#' odm_multi <- multimodal_ondisc_matrix(list(gene = odm_gene, grna = odm_grna))
#'
#' # delete gene_p_mito from cell covariate
#' odm_multi_no_gene_p_mito <- mutate_cell_covariates(odm_multi, gene_p_mito = NULL)
#' # add one column for cell covariate
#' odm_p_feature_expressed <- mutate_cell_covariates(odm_multi,
#' p_feature_expressed = gene_n_nonzero/(nrow(odm)[1]))
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

# covariate_ondisc_matrix functions
#' @rdname mutate-covariates
#' @export
setMethod("mutate_cell_covariates", signature("covariate_ondisc_matrix"), function(x, ...) {
  x@cell_covariates <- dplyr::mutate(x@cell_covariates, ...)
  return(x)
})

#' @rdname mutate-covariates
#' @export
setMethod("mutate_feature_covariates", signature("covariate_ondisc_matrix"), function(x, ...) {
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
#' @examples
#' # Install the `ondiscdata` package to run the examples.
#' # devtools::install_github("timothy-barry/ondiscdata")
#'
#' odm_gene_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_gene_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' odm_gene <- read_odm(odm_gene_fp, metadata_gene_fp)
#' odm_grna_fp <- system.file("extdata",
#' "odm/grna_assignment/matrix.odm", package = "ondiscdata")
#' metadata_grna_fp <- system.file("extdata",
#' "odm/grna_assignment/metadata.rds", package = "ondiscdata")
#' odm_grna <- read_odm(odm_grna_fp, metadata_grna_fp)
#'
#' odm_multi <- multimodal_ondisc_matrix(list(gene = odm_gene, grna = odm_grna))
#' odm_gene <- get_modality(odm_multi, "gene")
get_modality <- function(multimodal_mat, modality_name) {
  return(multimodal_mat@modalities[[modality_name]])
}


# 3. Get covariate matrices
#' Get covariate matrices.
#'
#' Functions to get cell-specific or feature-specific covariate matrices from  `covariate_ondisc_matrix` or `multimodal_ondisc_matrix` objects.
#' For `covariate_ondisc_matrix`, `get_cell_covariates` function can be used to get the cell-specific covariate matrix and `get_feature_covariates` can be used to get the feature-specific covariate matrix.
#' For`multimodal_ondisc_matrix`, `get_cell_covariates` function can be used to get the `global_cell_covariates` but `get_feature_covariates` is not supported.
#'
#' @name get-covariates
#' @param x an object of class `covariate_ondisc_matrix` or `multimodal_ondisc_matrix`.
#' @return a covariate matrix (in data frame form)
#' @examples
#' # Install the `ondiscdata` package to run the examples.
#' # devtools::install_github("timothy-barry/ondiscdata")
#'
#' ####################################
#' # EXAMPLE 1: covariate_ondisc_matrix
#' ####################################
#' odm_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' odm <- read_odm(odm_fp, metadata_fp)
#'
#' cell_covariates <- get_cell_covariates(odm)
#' feature_covariates <- get_feature_covariates(odm)
#'
#' #####################################
#' # EXAMPLE 2: multimodal_ondisc_matrix
#' #####################################
#' odm_gene_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_gene_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' odm_gene <- read_odm(odm_gene_fp, metadata_gene_fp)
#' odm_grna_fp <- system.file("extdata",
#' "odm/grna_assignment/matrix.odm", package = "ondiscdata")
#' metadata_grna_fp <- system.file("extdata",
#' "odm/grna_assignment/metadata.rds", package = "ondiscdata")
#' odm_grna <- read_odm(odm_grna_fp, metadata_grna_fp)
#' odm_multi <- multimodal_ondisc_matrix(list(gene = odm_gene, grna = odm_grna))
#'
#' cell_covariates_multi <- get_cell_covariates(odm_multi)
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

# covariate_ondisc_matrix functions
#' @rdname get-covariates
#' @export
setMethod("get_cell_covariates", signature("covariate_ondisc_matrix"), function(x) x@cell_covariates)

#' @rdname get-covariates
#' @export
setMethod("get_feature_covariates", signature("covariate_ondisc_matrix"), function(x) x@feature_covariates)

# multimodal_ondisc_matrix functions
#' @rdname get-covariates
#' @export
setMethod("get_cell_covariates", signature("multimodal_ondisc_matrix"), function(x) x@global_cell_covariates)

#' @rdname get-covariates
#' @export
setMethod("get_feature_covariates", signature("multimodal_ondisc_matrix"), function(x) {
  stop("get_feature_covariates not supported for objects of class multimodal_ondisc_matrix.")
})


# 4. Get `ondisc_matrix`
#'
#' Get the ondisc_matrix stored in a `covariate_ondisc_matrix.`
#'
#' @param covariate_odm a `covariate_ondisc_matrix` object
#'
#' @return the `ondisc_matrix` stored within the `covariate_ondisc_matrix` object
#' @export
#' @examples
#' # Install the `ondiscdata` package to run the examples.
#' # devtools::install_github("timothy-barry/ondiscdata")
#'
#' odm_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' cov_odm <- read_odm(odm_fp, metadata_fp)
#'
#' odm_mtx <- get_ondisc_matrix(cov_odm)
get_ondisc_matrix <- function(covariate_odm) {
  return(covariate_odm@ondisc_matrix)
}


####################################################
# ADD: select/extract cells, select/extract features
####################################################
