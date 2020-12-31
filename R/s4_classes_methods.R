# Class definition and methods for on_disc_matrix


# Classes and constructor
#########################
#' `on_disc_matrix` class constructor
#'
#' An on_disc_matrix represents a gene-by-cell expression matrix stored on disk. Use this function to obtain an `on_disc_matrix` from an already-initialized on-disk on_disc_matrix.h5 file.
#'
#' @slot h5_file path to an expression matrix (created by one of the initialization functions in this package) stored on disk
#' @slot cell_subset an integer vector storing the cells currently in use
#' @slot gene_subset an integer vector storing the genes currently in use
#'
#' @export on_disc_matrix
#' @examples
#' # NOTE: You must create the HDF5 file "example.h5" to run this example.
#' # Navigate to the help file of "create_on_disc_matrix_from_10x_mtx"
#' # (via ?create_on_disc_matrix_from_10x_mtx), and execute the code in the example.
#' odm_fp <- system.file("extdata", "example.h5", package = "ondisc")
#' if (odm_fp != "") { # if required file exists, ...
#' odm <- on_disc_matrix(h5_file = odm_fp)
#' }
on_disc_matrix <- setClass("on_disc_matrix",
                           representation(h5_file = "character", cell_subset = "integer", gene_subset = "integer"),
                           prototype(h5_file = NA_character_, cell_subset = NA_integer_, gene_subset = NA_integer_))

#' index_vec class
#' A union of the numeric, logical, and character classes. Used for indexing of on_disc_matrices. Idea borrowed from Matrix package.
setClassUnion("index_vec", members =  c("numeric", "logical", "character"))

# Basic information extraction methods
######################################

#' Return dimension
#'
#' Return dimension of on_disc_matrix.
#' @param x an on_disc_matrix
#' @export
#' @return an integer vector containing the dimension of the matrix
#' @examples
#' # NOTE: You must create the HDF5 file "example.h5" to run this example.
#' # Navigate to the help file of "create_on_disc_matrix_from_10x_mtx"
#' # (via ?create_on_disc_matrix_from_10x_mtx), and execute the code in the example.
#' odm_fp <- system.file("extdata", "example.h5", package = "ondisc")
#' if (odm_fp != "") { # if required file exists, ...
#' odm <- on_disc_matrix(h5_file = odm_fp)
#' dim(odm)
#' }
setMethod("dim", signature("on_disc_matrix"), function(x) get_dim(x))

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
#' odm <- on_disc_matrix(h5_file = odm_fp)
#' show(odm)
#' }
setMethod("show", signature = signature("on_disc_matrix"), function(object) {
  x_dim <- dim(object)
  cat(paste0("An on_disc_matrix with ", crayon::blue(x_dim[1]), " gene", if (x_dim[1] == 1) "" else "s" ," and ", crayon::blue(x_dim[2]), " cell", if (x_dim[2] == 1) "" else "s", ".\n"))
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
#' odm <- on_disc_matrix(h5_file = odm_fp)
#' head(odm)
#' }
setMethod("head", signature = signature("on_disc_matrix"), function(x) {
  x_dim <- dim(x)
  n_row_to_show <- min(5, x_dim[1])
  n_col_to_show <- min(6, x_dim[2])
  cat(paste0("Showing ", crayon::blue(n_row_to_show), " of ", crayon::blue(x_dim[1]), " gene", if (x_dim[1] == 1) "" else "s" ," and ", crayon::blue(n_col_to_show), " of ", crayon::blue(x_dim[2])," cell", if (n_col_to_show == 1) "" else "s", ":\n"))
  print(as.matrix(x[[1:n_row_to_show, 1:n_col_to_show]]))
})


# Subset methods and doc
########################

#' Subset an `on_disc_matrix` with `[`
#'
#' The user can pass logical, character, of numeric vectors to \code{`[`}. Character vectors correspond to gene IDs (for rows) and cell barcodes (for columns).
#' @param x An on_disc_matrix object
#' @param i Vector (numeric, logical, or character) indicating genes to keep
#' @param j Vector (numeric, logical, or character) indicating cells to keep
#' @param drop not used
#' @return A subset on_disc_matrix object.
#' @name subset-odm
#' @examples
#' # NOTE: You must create the HDF5 file "example.h5" to run this example.
#' # Navigate to the help file of "create_on_disc_matrix_from_10x_mtx"
#' # (via ?create_on_disc_matrix_from_10x_mtx), and execute the code in the example.
#' odm_fp <- system.file("extdata", "example.h5", package = "ondisc")
#' if (odm_fp != "") { # if required file exists, ...
#' odm <- on_disc_matrix(h5_file = odm_fp)
#' # keep cells 100-110
#' x <- odm[,100:110]
#' # keep genes ENSG00000260528, ENSG00000258908, and ENSG00000081913
#' x <- odm[c("ENSG00000260528", "ENSG00000258908", "ENSG00000081913"),]
#' # keep cells ACAGCCGCAGAAACCG and CTACTATAGTGTACCT
#' x <- odm[,c("ACAGCCGCAGAAACCG", "CTACTATAGTGTACCT")]
#' # keep all genes except ENSG00000167525 and ENSG00000235815
#' x <- odm[!(get_gene_ids(odm) %in% c("ENSG00000167525", "ENSG00000235815")),]
#' }
NULL

#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "on_disc_matrix", i = "missing", j = "missing", drop = "missing"),
          definition = function(x, i, j, drop) return(x))

#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "on_disc_matrix", i = "index_vec", j = "missing", drop = "missing"),
          definition = function(x, i, j) subset_by_gene_or_cell(x = x, idx = i, subset_on_cell = FALSE))

#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "on_disc_matrix", i = "missing", j = "index_vec", drop = "missing"),
          definition = function(x, i, j) subset_by_gene_or_cell(x = x, idx = j, subset_on_cell = TRUE))

#' @rdname subset-odm
#' @export
setMethod(f = "[",
          signature = signature(x = "on_disc_matrix", i = "index_vec", j = "index_vec", drop = "missing"),
          definition = function(x, i, j) {
            subset_by_gene_or_cell(x = x, idx = i, subset_on_cell = FALSE) %>% subset_by_gene_or_cell(x = ., idx = j, subset_on_cell = TRUE)
          })


# Extract expression data methods
#################################

#' Pull a submatrix into memory with `[[`
#'
#' The user can pass logical, character, of numeric vectors to \code{`[[`}. Character vectors correspond to gene IDs (for rows) and cell barcodes (for columns).
#' @param x An on_disc_matrix object
#' @param i Vector (numeric, logical, or character) indicating genes to keep
#' @param j Vector (numeric, logical, or character) indicating cells to keep
#' @return A matrix (of class Matrix)
#' @examples
#' # NOTE: You must create the HDF5 file "example.h5" to run this example.
#' # Navigate to the help file of "create_on_disc_matrix_from_10x_mtx"
#' # (via ?create_on_disc_matrix_from_10x_mtx), and execute the code in the example.
#' odm_fp <- system.file("extdata", "example.h5", package = "ondisc")
#' if (odm_fp != "") { # if required file exists, ...
#' odm <- on_disc_matrix(h5_file = odm_fp)
#' # keep cells 100-110
#' x <- odm[[,100:110]]
#' # keep genes ENSG00000260528, ENSG00000258908, and ENSG00000081913
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
          signature = signature(x = "on_disc_matrix", i = "missing", j = "missing"),
          definition = function(x, i, j) stop("Specify row or column indices to extract a sub-matrix."))

# 2. Extract by gene
#' @export
#' @rdname extract-odm
setMethod(f = "[[",
          signature = signature(x = "on_disc_matrix", i = "index_vec", j = "missing"),
          definition = function(x, i, j) subset_by_gene_or_cell(x = x, idx = i, subset_on_cell = FALSE) %>% extract_matrix())

# 3. Extract by cell
#' @export
#' @rdname extract-odm
setMethod(f = "[[",
          signature = signature(x = "on_disc_matrix", i = "missing", j = "index_vec"),
          definition = function(x, i, j) subset_by_gene_or_cell(x = x, idx = j, subset_on_cell = TRUE) %>% extract_matrix())

# 4. Extract by both gene and cell
#' @export
#' @rdname extract-odm
setMethod(f = "[[",
          signature = signature(x = "on_disc_matrix", i = "index_vec", j = "index_vec"),
          definition = function(x, i, j) subset_by_gene_or_cell(x = x, idx = i, subset_on_cell = FALSE) %>% subset_by_gene_or_cell(x = ., idx = j, subset_on_cell = TRUE) %>% extract_matrix())

# Summarize method, and more general apply/reduce methods
##########################################################
#' @export
setGeneric(name = "apply", def = function(X, MARGIN, FUN, ...) standardGeneric("apply"))

#' Apply an arbitrary function to all rows or columns
#' @param X an on_disc_matrix
#' @param MARGIN apply to rows (1) or columns (2)
#' @param FUN a function to apply
#' @param chunk_size number of rows or columns to load at a time; default 4000
#' @export
#' @examples
#' # NOTE: You must create the HDF5 file "example.h5" to run this example.
#' # Navigate to the help file of "create_on_disc_matrix_from_10x_mtx"
#' # (via ?create_on_disc_matrix_from_10x_mtx), and execute the code in the example.
#' odm_fp <- system.file("extdata", "example.h5", package = "ondisc")
#' if (odm_fp != "") { # if required file exists, ...
#' odm <- on_disc_matrix(h5_file = odm_fp)
#' gene_means_and_sds <- apply(X = odm,
#' MARGIN = 1,
#' FUN = function(r) c(mean = mean(r), sd = sd(r)),
#' chunk_size = 500)
#' }
setMethod(f = "apply",
          signature = signature("on_disc_matrix"),
          definition = function(X, MARGIN, FUN, chunk_size = 4000) {
            closure <- function(chunk) apply(X = as.matrix(chunk), MARGIN = MARGIN, FUN = FUN)
            on_disc_apply_across_chunks(x = X, col_apply = (MARGIN == 2), chunk_function = closure, chunk_size = chunk_size)
          })
