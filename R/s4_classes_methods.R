# Class definition and methods for ondisc_matrix


# Classes and constructor
#########################
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
#' # NOTE: You must create the HDF5 file "example.h5" to run this example.
#' # Navigate to the help file of "create_on_disc_matrix_from_10x_mtx"
#' # (via ?create_on_disc_matrix_from_10x_mtx), and execute the code in the example.
#' odm_fp <- system.file("extdata", "example.h5", package = "ondisc")
#' if (odm_fp != "") { # if required file exists, ...
#' odm <- ondisc_matrix(h5_file = odm_fp)
#' }
ondisc_matrix <- setClass("ondisc_matrix",
                           representation(h5_file = "character",
                                          logical_mat = "logical",
                                          cell_subset = "integer",
                                          cell_subset_order = "integer",
                                          feature_subset = "integer",
                                          feature_subset_order = "integer"),
                           prototype(h5_file = NA_character_,
                                     logical_mat = FALSE,
                                     cell_subset = NA_integer_,
                                     cell_subset_order = NA_integer_,
                                     feature_subset = NA_integer_,
                                     feature_subset_order = NA_integer_)
                          )

#' index_vec class
#' A union of the numeric, logical, and character classes. Used for indexing of on_disc_matrices. Idea borrowed from Matrix package.
# setClassUnion("index_vec", members =  c("numeric", "logical", "character"))

# Basic information extraction methods
######################################

#' Return dimension
#'
#' Return dimension of ondisc_matrix.
#' @param x an ondisc_matrix
#' @export
#' @return an integer vector containing the dimension of the matrix
#' @examples
#' # NOTE: You must create the HDF5 file "example.h5" to run this example.
#' # Navigate to the help file of "create_on_disc_matrix_from_10x_mtx"
#' # (via ?create_on_disc_matrix_from_10x_mtx), and execute the code in the example.
#' odm_fp <- system.file("extdata", "example.h5", package = "ondisc")
#' if (odm_fp != "") { # if required file exists, ...
#' odm <- ondisc_matrix(h5_file = odm_fp)
#' dim(odm)
#' }
setMethod("dim", signature("ondisc_matrix"), function(x) get_dim(x))

#' Print basic information to console
#'
#' Also display object class, number of rows, and number of columns.
#' @param object an on_dist_matrix to show
#' @return NULL
#' @export
#' @examples
#' # NOTE: You must create the HDF5 file "example.h5" to run this example.
#' # Navigate to the help file of "create_on_disc_matrix_from_10x_mtx"
#' # (via ?create_on_disc_matrix_from_10x_mtx), and execute the code in the example.
#' odm_fp <- system.file("extdata", "example.h5", package = "ondisc")
#' if (odm_fp != "") { # if required file exists, ...
#' odm <- ondisc_matrix(h5_file = odm_fp)
#' show(odm)
#' }
setMethod("show", signature = signature("ondisc_matrix"), function(object) {
  x_dim <- dim(object)
  cat(paste0("An ondisc_matrix with ", crayon::blue(x_dim[1]), " feature", if (x_dim[1] == 1) "" else "s" ," and ", crayon::blue(x_dim[2]), " cell", if (x_dim[2] == 1) "" else "s", ".\n"))
})

#' Print the first few rows and columns
#' @export
#' @param x an on_disc_mnatrix
#' @return NULL
#' @examples
#' # NOTE: You must create the HDF5 file "example.h5" to run this example.
#' # Navigate to the help file of "create_on_disc_matrix_from_10x_mtx"
#' # (via ?create_on_disc_matrix_from_10x_mtx), and execute the code in the example.
#' odm_fp <- system.file("extdata", "example.h5", package = "ondisc")
#' if (odm_fp != "") { # if required file exists, ...
#' odm <- ondisc_matrix(h5_file = odm_fp)
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
#' # NOTE: You must create the HDF5 file "example.h5" to run this example.
#' # Navigate to the help file of "create_on_disc_matrix_from_10x_mtx"
#' # (via ?create_on_disc_matrix_from_10x_mtx), and execute the code in the example.
#' odm_fp <- system.file("extdata", "example.h5", package = "ondisc")
#' if (odm_fp != "") { # if required file exists, ...
#' odm <- ondisc_matrix(h5_file = odm_fp)
#' # keep cells 100-110
#' x <- odm[,100:110]
#' # keep features ENSG00000260528, ENSG00000258908, and ENSG00000081913
#' x <- odm[c("ENSG00000260528", "ENSG00000258908", "ENSG00000081913"),]
#' # keep cells ACAGCCGCAGAAACCG and CTACTATAGTGTACCT
#' x <- odm[,c("ACAGCCGCAGAAACCG", "CTACTATAGTGTACCT")]
#' # keep all features except ENSG00000167525 and ENSG00000235815
#' x <- odm[!(get_feature_ids(odm) %in% c("ENSG00000167525", "ENSG00000235815")),]
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
#' # NOTE: You must create the HDF5 file "example.h5" to run this example.
#' # Navigate to the help file of "create_on_disc_matrix_from_10x_mtx"
#' # (via ?create_on_disc_matrix_from_10x_mtx), and execute the code in the example.
#' odm_fp <- system.file("extdata", "example.h5", package = "ondisc")
#' if (odm_fp != "") { # if required file exists, ...
#' odm <- ondisc_matrix(h5_file = odm_fp)
#' # keep cells 100-110
#' x <- odm[[,100:110]]
#' # keep features ENSG00000260528, ENSG00000258908, and ENSG00000081913
#' x <- odm[[c("ENSG00000260528", "ENSG00000258908", "ENSG00000081913"),]]
#' # keep cells ACAGCCGCAGAAACCG and CTACTATAGTGTACCT
#' x <- odm[[,c("ACAGCCGCAGAAACCG", "CTACTATAGTGTACCT")]]
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
