# Class definition and functions for on_disc_matrix


# Classes and constructor
#########################
#' The on_disc_matrix class
#'
#' An on_disc_matrix represents a gene-by-cell expression matrix stored on disk.
#'
#' @slot h5_file path to an expression matrix (created by one of the initialization functions in this package) stored on disk
#' @slot cell_subset an integer vector storing the cells currently in use
#' @slot gene_subset an integer vector storing the genes currently in use
#'
#' @export on_disc_matrix
on_disc_matrix <- setClass("on_disc_matrix",
                           representation(h5_file = "character", cell_subset = "integer", gene_subset = "integer"),
                           prototype(h5_file = NA_character_, cell_subset = NA_integer_, gene_subset = NA_integer_))

#' index_vec class
#' A union of the numeric, logical, and character classes. Used for indexing of on_disc_matrices. Idea borrowed from Matrix package.
setClassUnion("index_vec", members =  c("numeric", "logical", "character"))

# Basic information extraction methods
######################################

# 1. dim
#' dim
#' Print dimension of on_disc_matrix
#' @param x an on_disc_matrix
#' @export
#' @return an integer vector containing the dimension of the matrix
setMethod("dim", signature("on_disc_matrix"), function(x) get_dim(x))

# 2. show
#' print first few rows and columns of an on_disc_matrix to the console; also display object class
#' @param object an on_dist_matrix to show
#' @return NULL
#' @export
setMethod("show", signature = signature("on_disc_matrix"), function(object) {
  x_dim <- dim(object)
  cat(paste0("An on_disc_matrix with ", crayon::blue(x_dim[1]), " gene", if (x_dim[1] == 1) "" else "s" ," and ", crayon::blue(x_dim[2]), " cell", if (x_dim[2] == 1) "" else "s", ".\n"))
})

# 3. head
#' print first rows and columns of an on_disc_matrix
#' @export
#' @param x an on_disc_mnatrix
#' @return NULL
setMethod("head", signature = signature("on_disc_matrix"), function(x) {
  x_dim <- dim(x)
  n_row_to_show <- min(5, x_dim[1])
  n_col_to_show <- min(6, x_dim[2])
  cat(paste0("Showing ", crayon::blue(n_row_to_show), " of ", crayon::blue(x_dim[1]), " gene", if (x_dim[1] == 1) "" else "s" ," and ", crayon::blue(n_col_to_show), " of ", crayon::blue(x_dim[2])," cell", if (n_col_to_show == 1) "" else "s", ":\n"))
  print(as.matrix(x[[1:n_row_to_show, 1:n_col_to_show]]))
})


# Subset methods and doc
########################

#' Subset an on_disc_matrix
#' Subset an \linkS4class{on_disc_matrix} using the \code{`[`} operator.
#' @param x An \linkS4class{on_disc_matrix} object (Do not pass as an argument to \code{`[`}).
#' @param i Vector (numeric, logical, or character) indicating genes to keep.
#' @param j Vector (numeric, logical, or character) indicating cells to keep.
#' @param drop not used
#' @return A subset \linkS4class{on_disc_matrix} object.
#' @examples
#' # exp_mat_loc <- system.file("extdata", "on_disc_matrix_1.h5", package = "ondisc")
#' # if (exp_mat_loc != "") {
#' # x <- on_disc_matrix(h5_file = exp_mat_loc) # initialze an on_disc_matrix
#' # x_sub <- x[1:10,1:10] # subset first 10 genes and cells
#' # }
#' @name subset-odm
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

#' Extract an on_disc_matrix
#' Pull a submatrix of an \linkS4class{on_disc_matrix} into memory using the \code{`[[`} operator.
#' @param x An \linkS4class{on_disc_matrix} object (Do not pass as an argument to \code{`[[`}).
#' @param i Vector (numeric, logical, or character) indicating genes to keep.
#' @param j Vector (numeric, logical, or character) indicating cells to keep.
#' @return A matrix object (of class Matrix)
#' @examples
#' # exp_mat_loc <- system.file("extdata", "on_disc_matrix_1.h5", package = "ondisc")
#' # if (exp_mat_loc != "") {
#' # x <- on_disc_matrix(h5_file = exp_mat_loc) # initialze an on_disc_matrix
#' # x_mem <- x[[1:10,1:10]] # subset first 10 genes and cells
#' # }
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

# Col/Row sums methods, and more general apply/reduce methods
#############################################################
#' @export
setGeneric(name = "apply", def = function(X, MARGIN, FUN, ...) standardGeneric("apply"))

#' apply
#' apply a function to the rows or columns of an on_disc_matrix.
#' @param X an on_disc_matrix
#' @param MARGIN apply to rows (1) or columns (2)
#' @param FUN a function to apply
#' @param chunk_size number of rows or columns to load at a time; default 4000
#' @export
setMethod(f = "apply",
          signature = signature("on_disc_matrix"),
          definition = function(X, MARGIN, FUN, chunk_size = 4000) {
            closure <- function(chunk) apply(X = as.matrix(chunk), MARGIN = MARGIN, FUN = FUN)
            on_disc_apply_across_chunks(x = X, col_apply = (MARGIN == 2), chunk_function = closure, chunk_size = chunk_size)
          })
