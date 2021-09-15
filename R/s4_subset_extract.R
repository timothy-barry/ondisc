# Subset methods
################

#' Subset using the `[` operator.
#'
#' Apply the `[` operator to an `ondisc_matrix`, `covariate_ondisc_matrix`, or `multimodal_ondisc_matrix`
#' to subset the object. You can pass logical, character, or numeric vectors to `[`; character
#' vectors are assumed to refer to feature IDs (for rows) and cell barcodes (for columns).
#'
#' You can subset an `ondisc_matrix` and a `covariate_ondisc_matrix` by cell and/or feature. You can subset a
#' `multimodal_ondisc_matrix` by cell only (because the features differ across modalities).
#'
#' @param x an `ondisc_matrix`, `covariate_ondisc_matrix`, or `multimodal_ondisc_matrix` object.
#' @param i a vector (numeric, logical, or character) indicating features to keep.
#' @param j a vector (numeric, logical, or character) indicating cells to keep.
#' @param drop not used
#' @return An appropriately subset object of the same class as `x`.
#' @name subset-odm
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
#'
#' # numeric vector
#' odm_numeric <- odm[1:10, 1:10]
#' # logical vector
#' set.seed(1)
#' nrow_logical <- sample(x = c(TRUE, FALSE), size = nrow(odm), replace = TRUE)
#' ncol_logical <- sample(x = c(TRUE, FALSE), size = ncol(odm), replace = TRUE)
#' odm_logical <- odm[nrow_logical, ncol_logical]
#' # character vectors for feature IDs and cell barcodes
#' feature_ids_subset <- get_feature_ids(odm)[1:10]
#' cell_barcodes_subset <- get_cell_barcodes(odm)[1:10]
#' odm_char <- odm[feature_ids_subset,cell_barcodes_subset]
#'
#' #####################################
#' # EXAMPLE 2: multimodal_ondisc_matrix
#' #####################################
#' odm_gene_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_gene_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' odm_gene <- read_odm(odm_gene_fp, metadata_gene_fp)
#' odm_gRNA_fp <- system.file("extdata", "odm/gRNA/matrix.odm", package = "ondiscdata")
#' metadata_gRNA_fp <- system.file("extdata", "odm/gRNA/metadata.rds", package = "ondiscdata")
#' odm_gRNA <- read_odm(odm_gRNA_fp, metadata_gRNA_fp)
#' odm_multi <- multimodal_ondisc_matrix(list(gene = odm_gene, gRNA = odm_gRNA))
#'
#' # numeric vector
#' odm_multi_numeric <- odm_multi[, 1:10]
#' # logical vector
#' set.seed(1)
#' ncol_logical <- sample(x = c(TRUE, FALSE), size = ncol(odm_multi), replace = TRUE)
#' odm_multi_logical <- odm_multi[, ncol_logical]
#' # character vectors for cell barcodes
#' cell_barcodes_subset <- get_cell_barcodes(odm_multi)[1:10]
#' odm_multi_cell_barcodes <- odm[,cell_barcodes_subset]
NULL

# 1. subset nothing, odm
#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "ondisc_matrix", i = "missing", j = "missing", drop = "missing"),
          definition = function(x, i, j, drop) return(x))

# 2. subset cell, odm
#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "ondisc_matrix", i = "ANY", j = "missing", drop = "missing"),
          definition = function(x, i, j) subset_by_feature_or_cell(x = x, idx = i, subset_on_cell = FALSE, subset_string_arrays = TRUE))

# 3. subset feature, odm
#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "ondisc_matrix", i = "missing", j = "ANY", drop = "missing"),
          definition = function(x, i, j) subset_by_feature_or_cell(x = x, idx = j, subset_on_cell = TRUE, subset_string_arrays = TRUE))

# 4. subset both cell and feature, odm
#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "ondisc_matrix", i = "ANY", j = "ANY", drop = "missing"),
          definition = function(x, i, j) {
            subset_by_feature_or_cell(x = x, idx = i, subset_on_cell = FALSE, subset_string_arrays = TRUE) %>% subset_by_feature_or_cell(x = ., idx = j, subset_on_cell = TRUE, subset_string_arrays = TRUE)
          })

# 5. subset both cell and feature, m_odm
#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "covariate_ondisc_matrix", i = "ANY", j = "ANY", drop = "missing"),
          definition = function(x, i, j, drop) {
            x@ondisc_matrix <- x@ondisc_matrix[i,j]
            x@cell_covariates <- x@cell_covariates[j,,drop=FALSE]
            x@feature_covariates <- x@feature_covariates[i,,drop=FALSE]
            return(x)
          })

# 6. subset cell, m_odm
#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "covariate_ondisc_matrix", i = "ANY", j = "missing", drop = "missing"),
          definition = function(x, i, j, drop) {
            x@ondisc_matrix <- x@ondisc_matrix[i,]
            x@feature_covariates <- x@feature_covariates[i,,drop=FALSE]
            return(x)
          })

# 7. subset feature, m_odm
#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "covariate_ondisc_matrix", i = "missing", j = "ANY", drop = "missing"),
          definition = function(x, i, j, drop) {
            x@ondisc_matrix <- x@ondisc_matrix[,j]
            x@cell_covariates <- x@cell_covariates[j,,drop=FALSE]
            return(x)
          })

# 8. subset nothing, m_odm
#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "covariate_ondisc_matrix", i = "missing", j = "missing", drop = "missing"),
          definition = function(x, i, j, drop) return(x))

# 9. subset nothing, mm_odm
#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "multimodal_ondisc_matrix", i = "missing", j = "missing", drop = "missing"),
          definition = function(x, i, j, drop) return(x))

# 10. subset feature, mm_odm
#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "multimodal_ondisc_matrix", i = "missing", j = "ANY", drop = "missing"),
          definition = function(x, i, j, drop) {
            n_modalities <- length(x@modalities)
            for (ctr in seq(1, n_modalities)) x@modalities[[ctr]] <- x@modalities[[ctr]][,j]
            x@global_cell_covariates <- x@global_cell_covariates[j,]
            return(x)
          })

# 11. subset cell, mm_odm
#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "multimodal_ondisc_matrix", i = "ANY"),
          definition = function(x, i, j, drop) stop(
            "The operation [i,] is invalid, because multimodal_ondisc_matrix objects can
            be subset only by cell, not feature. Use [,i] to subset by cell instead."))


#################
# Extract methods
#################

#' Pull a submatrix into memory using the `[[` operator.
#'
#' Apply the `[[` operator to an `ondisc_matrix` or `covariate_ondisc_matrix` to pull a submatrix into memory. You can pass logical, character,
#' or numeric vectors to `[[`; character vectors are assumed to refer to feature IDs (for rows) and cell barcodes
#' (for columns).
#'
#' One can remember the difference between `[` and `[[` by recalling R lists: `[` is used to subset a list, and
#' `[[` is used to access elements stored *inside* a list. Similarly, `[` is used to subset an `ondisc_matrix`, and
#' `[[` is used to access a submatrix usable within R.
#'
#' @param x an `ondisc_matrix` object.
#' @param i a vector (numeric, logical, or character) indicating features to pull into memory.
#' @param j a vector (numeric, logical, or character) indicating cells to pull into memory.
#' @return a matrix (as implemented by the Matrix package).
#' @name extract-odm
#' @examples
#' # Use `ondiscdata` package for the examlpes, please install the package before running the examples
#' # install.packages("devtools")
#' # devtools::install_github("Katsevich-Lab/ondiscdata")
#'
#' odm_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' odm <- read_odm(odm_fp, metadata_fp)
#'
#' # pull a single feature
#' single_feature <- odm[["ENSG00000180098.9",]]
#' # pull a single cell
#' single_cell <- odm[[,"GGAGCAAGTACTTAGC_2_2"]]
#' # numeric vector
#' mtx_numeric <- odm[[1:10, 1:10]]
#' # logical vector
#' set.seed(1)
#' nrow_logical <- sample(x = c(TRUE, FALSE), size = nrow(odm), replace = TRUE)
#' ncol_logical <- sample(x = c(TRUE, FALSE), size = ncol(odm), replace = TRUE)
#' mtx_logical <- odm[[nrow_logical, ncol_logical]]
#' # character vectors for feature IDs and cell barcodes
#' feature_ids_subset <- get_feature_ids(odm)[1:10]
#' cell_barcodes_subset <- get_cell_barcodes(odm)[1:10]
#' mtx_char <- odm[[feature_ids_subset, cell_barcodes_subset]]
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

# 5. metadata_odm
#' @export
#' @rdname extract-odm
setMethod(f = "[[",
          signature = signature(x = "covariate_ondisc_matrix", i = "ANY", j = "ANY"),
          definition = function(x, i, j)
          x@ondisc_matrix[[i, j]]
)

#' @export
#' @rdname extract-odm
setMethod(f = "[[",
          signature = signature(x = "covariate_ondisc_matrix", i = "missing", j = "ANY"),
          definition = function(x, i, j)
            x@ondisc_matrix[[, j]]
)

#' @export
#' @rdname extract-odm
setMethod(f = "[[",
          signature = signature(x = "covariate_ondisc_matrix", i = "ANY", j = "missing"),
          definition = function(x, i, j)
            x@ondisc_matrix[[i, ]]
)

#' @export
#' @rdname extract-odm
setMethod(f = "[[",
          signature = signature(x = "covariate_ondisc_matrix", i = "missing", j = "missing"),
          definition = function(x, i, j)
            x@ondisc_matrix[[, ]]
)


# 6. multimodal_odm
#' @export
#' @rdname extract-odm
setMethod(f = "[[",
          signature = signature(x = "multimodal_ondisc_matrix"),
          definition = function(x, i, j) stop(
            "You cannot use the [[,]] operator on this object because it is a multimodal_ondisc_matrix,
            not an ondisc_matrix. To access an ondisc_matrix stored within this object, (i) access the
            `modalities` field using @, (ii) select a modality using $, and (iii) access the `ondisc_matrix`
            field using @."
          ))
