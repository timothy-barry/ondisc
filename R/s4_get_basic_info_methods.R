# Dimension
###########

#' Get dimension
#'
#' Return the dimension of an `ondisc_matrix`, `covariate_ondisc_matrix`, or `multimodal_ondisc_matrix`.
#'
#' @param x an `ondisc_matrix`, `covariate_ondisc_matrix`, or `multimodal_ondisc_matrix`.
#' @name dim
#' @return If `x` is an `ondisc_matrix` or `covariate_ondisc_matrix`, length-two integer vector containing
#' the dimension of `x`; if `x` is a `multimodal_ondisc_matrix`, a list of integer vectors containing the dimensions
#' of the constituent modalities of `x`.
#'
#' @examples
#' # Use `ondiscdata` package for the examlpes, please install the package before running the examples
#' # install.packages("devtools")
#' # devtools::install_github("Katsevich-Lab/ondiscdata")
#'
#' ####################################
#' # EXAMPLE 1: covariate_ondisc_matrix
#' ####################################
#' odm_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' odm <- read_odm(odm_fp, metadata_fp)
#' dim(odm)
#' nrow(odm)
#' ncol(odm)
#'
#' #####################################
#' # EXAMPLE 2: multimodal_ondisc_matrix
#' #####################################
#' # Load odm from package
#' odm_gene_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_gene_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' odm_gene <- read_odm(odm_gene_fp, metadata_gene_fp)
#' odm_gRNA_fp <- system.file("extdata", "odm/gRNA/matrix.odm", package = "ondiscdata")
#' metadata_gRNA_fp <- system.file("extdata", "odm/gRNA/metadata.rds", package = "ondiscdata")
#' odm_gRNA <- read_odm(odm_gRNA_fp, metadata_gRNA_fp)
#'
#' odm_multi <- multimodal_ondisc_matrix(list(gene = odm_gene, gRNA = odm_gRNA))
#' dim(odm_multi)
#' nrow(odm_multi)
#' ncol(odm_multi)
NULL

#' @export
#' @rdname dim
setMethod("dim", signature("ondisc_matrix"), function(x) get_dim(x))

#' @export
#' @rdname dim
setMethod("dim", signature("covariate_ondisc_matrix"), function(x) get_dim(x@ondisc_matrix))

#' @export
#' @rdname dim
setMethod("dim", signature("multimodal_ondisc_matrix"), function(x) lapply(x@modalities, dim))

#' @export
#' @rdname dim
setGeneric("ncol", function(x) standardGeneric("ncol"))

#' @export
#' @rdname dim
setGeneric("nrow", function(x) standardGeneric("nrow"))

#' @export
#' @rdname dim
setMethod("ncol", signature("multimodal_ondisc_matrix"), function(x) ncol(x@modalities[[1]]))

#' @export
#' @rdname dim
setMethod("nrow", signature("multimodal_ondisc_matrix"), function(x) lapply(x@modalities, nrow))


# Show
######

#' Print basic information to the console
#'
#' @param object an object of class ondisc_matrix, covaraite_ondisc_matrix, or multimodal_ondisc_matrix
#' @name show
#' @return NULL; called for printing
NULL

#' @rdname show
#' @export
setMethod("show", signature = signature("ondisc_matrix"), function(object) {
  x_dim <- dim(object)
  cat(paste0("An ondisc_matrix with ", crayon::blue(x_dim[1]), " feature", if (x_dim[1] == 1) "" else "s" ," and ", crayon::blue(x_dim[2]), " cell", if (x_dim[2] == 1) "" else "s", ".\n"))
})

#' @rdname show
#' @export
setMethod("show", signature = signature("covariate_ondisc_matrix"), function(object) {
  cell_covariates <- colnames(object@cell_covariates)
  feature_covariates <- colnames(object@feature_covariates)
  cat("A covariate_ondisc_matrix with the following components:\n")
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


# Head
######

#' Head
#'
#' Print the first few rows and columns of an `ondisc_matrix`, `covariate_ondisc_matrix`, or `multimodal_ondisc_matrix`.
#'
#' @param x an `ondisc_matrix`, `covariate_ondisc_matrix`, or `multimodal_ondisc_matrix`.
#' @name head
#' @return NULL; called for printing
#' @examples
#' # Use `ondiscdata` package for the examlpes, please install the package before running the examples
#' # install.packages("devtools")
#' # devtools::install_github("Katsevich-Lab/ondiscdata")
#'
#' ####################################
#' # EXAMPLE 1: covariate_ondisc_matrix
#' ####################################
#' # Load odm from package
#' odm_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' odm <- read_odm(odm_fp, metadata_fp)
#' head(odm)
#'
#' #####################################
#' # EXAMPLE 2: multimodal_ondisc_matrix
#' #####################################
#' # Load odm from package
#' odm_gene_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_gene_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' odm_gene <- read_odm(odm_gene_fp, metadata_gene_fp)
#' odm_gRNA_fp <- system.file("extdata", "odm/gRNA/matrix.odm", package = "ondiscdata")
#' metadata_gRNA_fp <- system.file("extdata", "odm/gRNA/metadata.rds", package = "ondiscdata")
#' odm_gRNA <- read_odm(odm_gRNA_fp, metadata_gRNA_fp)
#'
#' odm_multi <- multimodal_ondisc_matrix(list(gene = odm_gene, gRNA = odm_gRNA))
#' head(odm_multi)
NULL

#' @export
#' @rdname head
setMethod("head", signature = signature("ondisc_matrix"), function(x) {
  x_dim <- dim(x)
  n_row_to_show <- min(5, x_dim[1])
  n_col_to_show <- min(6, x_dim[2])
  cat(paste0("Showing ", crayon::blue(n_row_to_show), " of ", crayon::blue(x_dim[1]), " features", if (x_dim[1] == 1) "" else "s" ," and ", crayon::blue(n_col_to_show), " of ", crayon::blue(x_dim[2])," cell", if (n_col_to_show == 1) "" else "s", ":\n"))
  print(as.matrix(x[[1:n_row_to_show, 1:n_col_to_show]]))
})

#' @export
#' @rdname head
setMethod("head", signature = signature("covariate_ondisc_matrix"), function(x) {
  head(x@ondisc_matrix)
})

#' @export
#' @rdname head
setMethod("head", signature = signature("multimodal_ondisc_matrix"), function(x) {
  for (i in seq(1, length(x@modalities))) {
    paste0(i, ". ", crayon::blue(names(x@modalities[i])), " ") %>% cat()
    head(x@modalities[[i]])
  }
})


# Get feature names, feature ids, and cell barcodes
###################################################

#' Get cell barcodes, feature names, and feature IDs
#'
#' Obtain cell barcodes, feature names, and feature IDs of an `ondisc_matrix`, `covariate_ondisc_matrix`,
#' or `multimodal_ondisc_matrix`.
#'
#' The following functions can be used to obtain feature and cell identifiers:
#' - `get_cell_barcodes`: return the cell barcodes.
#' - `get_feature_names`: return the feature names.
#' - `get_feature_ids`: return the IDs of the features.
#'
#'
#' In general, these functions return a character vector containing the requested identifiers. When
#' `get_feature_names` or `get_feature_ids` is called on a `multimodal_ondisc_matrix`, the function instead
#' returns a list containing the feature names and feature IDs, respectively, of the modalities contained
#' within the `multimodal_ondisc_matrix`.
#'
#' @name get-names
#' @param x an object of class `ondisc_matrix`, `covariate_ondisc_matrix`, or `multimodal_ondisc_matrix`.
#' @return A character vector or list of character vectors containing the requested identifiers.
#'
#' @examples
#' # Use `ondiscdata` package for the examlpes, please install the package before running the examples
#' # install.packages("devtools")
#' # devtools::install_github("Katsevich-Lab/ondiscdata")
#'
#' ####################################
#' # EXAMPLE 1: covariate_ondisc_matrix
#' ####################################
#' # Load odm from package
#' odm_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' odm <- read_odm(odm_fp, metadata_fp)
#' feature_ids <- get_feature_ids(odm)
#' feature_names <- get_feature_names(odm)
#' cell_barcodes <- get_cell_barcodes(odm)
#'
#' #####################################
#' # EXAMPLE 2: multimodal_ondisc_matrix
#' #####################################
#' # Load odm from package
#' odm_gene_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_gene_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' odm_gene <- read_odm(odm_gene_fp, metadata_gene_fp)
#' odm_gRNA_fp <- system.file("extdata", "odm/gRNA/matrix.odm", package = "ondiscdata")
#' metadata_gRNA_fp <- system.file("extdata", "odm/gRNA/metadata.rds", package = "ondiscdata")
#' odm_gRNA <- read_odm(odm_gRNA_fp, metadata_gRNA_fp)
#'
#' odm_multi <- multimodal_ondisc_matrix(list(gene = odm_gene, gRNA = odm_gRNA))
#' feature_ids <- get_feature_ids(odm_multi)
#' feature_names <- get_feature_names(odm_multi)
#' cell_barcodes <- get_cell_barcodes(odm_multi)
NULL

#' @export
#' @rdname get-names
setGeneric("get_feature_ids", function(x) standardGeneric("get_feature_ids"))
#' @export
#' @rdname get-names
setGeneric("get_feature_names", function(x) standardGeneric("get_feature_names"))
#' @export
#' @rdname get-names
setGeneric("get_cell_barcodes", function(x) standardGeneric("get_cell_barcodes"))

#' @rdname get-names
#' @export
setMethod("get_feature_ids", signature("ondisc_matrix"), function(x) get_names(x, "feature_ids"))
#' @rdname get-names
#' @export
setMethod("get_feature_names", signature("ondisc_matrix"), function(x) get_names(x, "feature_names"))
#' @rdname get-names
#' @export
setMethod("get_cell_barcodes", signature("ondisc_matrix"), function(x) get_names(x, "cell_barcodes"))

#' @rdname get-names
#' @export
setMethod("get_feature_ids", signature("covariate_ondisc_matrix"), function(x) get_feature_ids(x@ondisc_matrix))
#' @rdname get-names
#' @export
setMethod("get_feature_names", signature("covariate_ondisc_matrix"), function(x) get_feature_names(x@ondisc_matrix))
#' @rdname get-names
#' @export
setMethod("get_cell_barcodes", signature("covariate_ondisc_matrix"), function(x) get_cell_barcodes(x@ondisc_matrix))

#' @rdname get-names
#' @export
setMethod("get_feature_ids", signature("multimodal_ondisc_matrix"), function(x) lapply(x@modalities, get_feature_ids))
#' @rdname get-names
#' @export
setMethod("get_feature_names", signature("multimodal_ondisc_matrix"), function(x) lapply(x@modalities, get_feature_names))
#' @rdname get-names
#' @export
setMethod("get_cell_barcodes", signature("multimodal_ondisc_matrix"), function(x) get_cell_barcodes(x@modalities[[1]]))
