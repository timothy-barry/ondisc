# Subset methods
################

#' Subset using the `[` operator.
#'
#' Apply the `[` operator to an `ondisc_matrix`, `metadata_ondisc_matrix`, or `multimodal_ondisc_matrix`
#' to subset the object. You can pass logical, character, or numeric vectors to `[`; character
#' vectors are assumed to refer to feature IDs (for rows) and cell barcodes (for columns).
#'
#' You can subset an `ondisc_matrix` and a `metadata_ondisc_matrix` by cell and/or feature. You can subset a
#' `multimodal_ondisc_matrix` by cell only (because the features differ across modalities).
#'
#' @param x an `ondisc_matrix`, `metadata_ondisc_matrix`, or `multimodal_ondisc_matrix` object.
#' @param i a vector (numeric, logical, or character) indicating features to keep.
#' @param j a vector (numeric, logical, or character) indicating cells to keep.
#' @param drop not used
#' @return An appropriately subset object of the same class as `x`.
#' @name subset-odm
#' @examples
#' # NOTE: You must create the RDS files "expressions.rds" and
#' # "perturbations.rds" to run this example. Navigate to the help file of
#' # "create_ondisc_matrix_from_mtx" (via ?create_ondisc_matrix_from_mtx),
#' # and execute both code blocks.
#'
#' # subset an ondisc_matrix
#' h5_fp <- paste0(tempdir(), "/expressions.h5")
#' if (file.exists(h5_fp)) {
#' odm <- ondisc_matrix(h5_file = h5_fp)
#' # keep cells 100-110
#' x <- odm[,100:110]
#' # keep all cells except 50, 100, 150
#' x <- odm[,-c(50, 100, 150)]
#' # keep genes ENSG00000188305, ENSG00000257284, and ENSG00000251655:
#' x <- odm[c("ENSG00000188305", "ENSG00000257284", "ENSG00000251655"),]
#' # keep the cells CTTAGGACACTGGCGT-1 and AAAGGATTCACATCAG-1:
#' x <- odm[,c("CTTAGGACACTGGCGT-1", "AAAGGATTCACATCAG-1")]
#' # keep all genes except ENSG00000188305 and ENSG00000257284
#' x <- odm[!(get_feature_ids(odm) %in% c("ENSG00000188305", "ENSG00000257284")),]
#' }
#'
#' # subset a metadata_ondic_matrix
#' expressions_fp <- paste0(tempdir(), "/expressions.rds")
#' if (file.exists(expressions_fp)) {
#' expressions <- readRDS(expressions_fp)
#' # keep cells 100-110
#' x <- expressions[,100:110]
#' # keep genes ENSG00000188305, ENSG00000257284, and ENSG00000251655
#' x <- expressions[c("ENSG00000188305", "ENSG00000257284", "ENSG00000251655"),]
#' }
#'
#' # subset a multimodal ondisc_matrix
#' expression_fp <- paste0(tempdir(), "/expressions.rds")
#' perturbations_fp <- paste0(tempdir(), "/perturbations.rds")
#' if (file.exists(expression_fp) && file.exists(perturbations_fp)) {
#'     expressions <- readRDS(expression_fp)
#'     perturbations <- readRDS(expression_fp)
#'     crispr_experiment <- multimodal_ondisc_matrix(list(expressions = expressions,
#'     perturbations = perturbations))
#'     # Keep all cells except 10,100, and 105.
#'     x <- crispr_experiment[,-c(10,100,105)]
#'     # Keep the first 5 cells
#'     x <- crispr_experiment[,1:5]
#' }
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
          definition = function(x, i, j) subset_by_feature_or_cell(x = x, idx = i, subset_on_cell = FALSE))

# 3. subset feature, odm
#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "ondisc_matrix", i = "missing", j = "ANY", drop = "missing"),
          definition = function(x, i, j) subset_by_feature_or_cell(x = x, idx = j, subset_on_cell = TRUE))

# 4. subset both cell and feature, odm
#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "ondisc_matrix", i = "ANY", j = "ANY", drop = "missing"),
          definition = function(x, i, j) {
            subset_by_feature_or_cell(x = x, idx = i, subset_on_cell = FALSE) %>% subset_by_feature_or_cell(x = ., idx = j, subset_on_cell = TRUE)
          })

# 5. subst both cell and feature, m_odm
#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "metadata_ondisc_matrix", i = "ANY", j = "ANY", drop = "missing"),
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
          signature = signature(x = "metadata_ondisc_matrix", i = "ANY", j = "missing", drop = "missing"),
          definition = function(x, i, j, drop) {
            x@ondisc_matrix <- x@ondisc_matrix[i,]
            x@feature_covariates <- x@feature_covariates[i,,drop=FALSE]
            return(x)
          })

# 7. subset feature, m_odm
#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "metadata_ondisc_matrix", i = "missing", j = "ANY", drop = "missing"),
          definition = function(x, i, j, drop) {
            x@ondisc_matrix <- x@ondisc_matrix[,j]
            x@cell_covariates <- x@cell_covariates[j,,drop=FALSE]
            return(x)
          })

# 8. subset nothing, m_odm
#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "metadata_ondisc_matrix", i = "missing", j = "missing", drop = "missing"),
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
#' Apply the `[[` operator to an `ondisc_matrix` to pull a submatrix into memory. You can pass logical, character,
#' or numeric vectors to `[[`; character vectors are assumed to refer to feature IDs (for rows) and cell barcodes
#' (for columns).
#'
#' You can apply `[[` to `ondisc_matrix` objects only. You cannot apply `[[` to `metadata_ondisc_matrix` or
#' `multimodal_ondisc_matrix` objects, because in the latter case the data to be accessed is ambiguous.
#'
#' You can remember the difference between `[` and `[[` by thinking about R lists: `[` is used to subset a list, and
#' `[[` is used to access elements stored *inside* a list. Similarly, `[` is used to subset an `ondisc_matrix`, and
#' `[[` is used to access a submatrix usable within R.
#'
#' @param x an `ondisc_matrix` object.
#' @param i a vector (numeric, logical, or character) indicating features to pull into memory.
#' @param j a vector (numeric, logical, or character) indicating cells to pull into memory.
#' @return a matrix (as implemented by the Matrix package).
#' @examples
#' # NOTE: You must create the HDF5 file "expressions.h5" to run this example.
#' # Navigate to the help file of "create_ondisc_matrix_from_mtx"
#' # (via ?create_ondisc_matrix_from_mtx), and execute the code in the first code block.
#'
#' h5_fp <- paste0(tempdir(), "/expressions.h5")
#' if (file.exists(h5_fp)) {
#' odm <- ondisc_matrix(h5_file = h5_fp)
#' # extract cells 100-110:
#' x <- odm[[,100:110]]
#' # extract genes ENSG00000188305, ENSG00000257284, ENSG00000251655:
#' x <- odm[[c("ENSG00000188305", "ENSG00000257284", "ENSG00000251655"),]]
#' # extract cells CTTAGGACACTGGCGT-1 and AAAGGATTCACATCAG-1:
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

# 5. metadata_odm
#' @export
#' @rdname extract-odm
setMethod(f = "[[",
          signature = signature(x = "metadata_ondisc_matrix"),
          definition = function(x, i, j) stop(
            "You cannot use the [[,]] operator on this object because it is a metadata_ondisc_matrix,
            not an ondisc_matrix. To access the ondisc_matrix stored within this object, use the @
            symbol and access the field `ondisc_matrix.`")
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
