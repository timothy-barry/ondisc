# odm class
setClass("odm",
         slots = list(h5_file = "character",
                      dimension = "integer",
                      feature_ids = "character",
                      ptr = "externalptr",
                      integer_id = "integer"))


#' Initialize an `odm` object
#'
#' `initialize_odm_from_backing_file()` initializes an `odm` object from a backing `.odm` file.
#'
#' `initialize_odm_from_backing_file()` is portable: the user can create an `.odm` file (via `create_odm_from_cellranger()` or `create_odm_from_r_matrix()`) on one computer, transfer the `.odm` file to another computer, and then load the `.odm` file (via `initialize_odm_from_backing_file()`) on the second computer.
#'
#' @param odm_file file path to a backing `.odm` file.
#'
#' @return an `odm` object
#' @export
#' @examples
#' library(sceptredata)
#' directories_to_load <- paste0(
#'  system.file("extdata", package = "sceptredata"),
#'  "/highmoi_example/gem_group_", c(1, 2)
#' )
#' directory_to_write <- tempdir()
#' out_list <- create_odm_from_cellranger(
#'   directories_to_load = directories_to_load,
#'   directory_to_write = directory_to_write,
#' )
#' gene_odm <- out_list$gene
#' gene_odm
#'
#' # delete the gene_odm object
#' rm(gene_odm)
#'
#' # reinitialize the gene_odm object
#' gene_odm <- initialize_odm_from_backing_file(
#'   paste0(tempdir(), "/gene.odm")
#' )
#' gene_odm
initialize_odm_from_backing_file <- function(odm_file) {
  # 0. checks on odm_file
  odm_file <- expand_tilde(odm_file)
  if (!file.exists(odm_file)) {
    stop(odm_file, " does not exist.")
  }
  if (!(methods::is(odm_file, "character") && length(odm_file) == 1)) {
    stop(odm_file, " must be a single string.")
  }

  # initialize the odm
  out <- methods::new("odm")
  h5_file <- odm_file
  dimension <- read_integer_vector(h5_file, "dimension", 2L)
  integer_id <- read_integer_vector(h5_file, "integer_id", 1L)
  feature_ids <- read_feature_ids(h5_file, dimension[1])
  ptr <- read_row_ptr(h5_file, dimension[1])
  out@h5_file <- h5_file
  out@dimension <- dimension
  out@integer_id <- integer_id
  out@feature_ids <- feature_ids
  out@ptr <- ptr
  return(out)
}


#' Load a row of an `odm` object into memory
#'
#' The operator `[` loads a specified row of an `odm` object into memory. The `odm` object can be indexed either by integer or feature ID.
#'
#' @param x an object of class `odm`
#' @param i the index of the row to load into memory. `i` can be either an integer (specifying the integer-based index of the row to load into memory) or a string (specifying the feature ID of the row to load into memory).
#' @param j not used
#' @param drop not used
#' @export
#' @return an expression vector (of class `numeric`)
#' @examples
#' library(sceptredata)
#' directories_to_load <- paste0(
#'  system.file("extdata", package = "sceptredata"),
#'  "/highmoi_example/gem_group_", c(1, 2)
#' )
#' directory_to_write <- tempdir()
#' out_list <- create_odm_from_cellranger(
#'   directories_to_load = directories_to_load,
#'   directory_to_write = directory_to_write,
#' )
#' gene_odm <- out_list$gene
#' # extract rows into memory by index and ID
#' v1 <- gene_odm[10L,]
#' v2 <- gene_odm["ENSG00000173825",]
setMethod(f = "[",
          signature = signature(x = "odm", i = "ANY", j = "missing", drop = "missing"),
          definition = function(x, i, j, drop) {
            if (length(i) != 1L) stop("The index must be a single element (i.e., not a vector).")
            if (methods::is(i, "numeric")) {
              n_features <- x@dimension[1]
              if (i <= 0 || i > n_features) {
                stop("The index must be an integer in the range [1, ", n_features, "].")
              }
            } else if (methods::is(i, "character")) {
              idx <- which(i == x@feature_ids)
              if (length(idx) != 1) {
                stop("The feature '", i, "' is not contained within the odm.")
              }
              i <- idx
            } else {
              stop("The row index must be of type numeric or character.")
            }
            load_row_cpp(file_name_in = x@h5_file, f_row_ptr = x@ptr,
                         row_idx = i - 1L, n_cells = x@dimension[2])
          })


setMethod(f = "show", signature = "odm", definition = function(object) {
  cat(paste0("An object of class ", crayon::blue("odm"), " with the following attributes:",
             "\n\t\U2022 ", crayon::blue(object@dimension[1]), " features",
             "\n\t\U2022 ", crayon::blue(object@dimension[2]), " cells",
             "\n\t\U2022 Backing file: ", crayon::blue(object@h5_file)))
})


#' Return the rownames of an `odm` object
#'
#' `rownames()` returns the rownames (i.e., feature IDs) of an `odm` object.
#'
#' @param x an object of class `odm`
#'
#' @return the rownames of an `odm` object
#' @aliases rownames
#' @rdname rownames-odm-method
#' @export
#' @examples
#' library(sceptredata)
#' directories_to_load <- paste0(
#'  system.file("extdata", package = "sceptredata"),
#'  "/highmoi_example/gem_group_", c(1, 2)
#' )
#' directory_to_write <- tempdir()
#' out_list <- create_odm_from_cellranger(
#'   directories_to_load = directories_to_load,
#'   directory_to_write = directory_to_write,
#' )
#' gene_odm <- out_list$gene
#' # return the rownames
#' rownames(gene_odm) |> head()
setMethod(f = "dimnames", signature = "odm", definition = function(x) list(x@feature_ids, NULL))


#' Return the number of columns and rows of an `odm` object
#'
#'  `ncol()` and `nrow()` return the number of rows (i.e., features) and columns (i.e., cells), respectively, contained within an `odm` object. `dim()` returns an integer vector of length two whose first and second entries, respectively, indicate the number of rows and columns in an `odm` object.
#'
#' @param x an object of class `odm`
#'
#' @return the dimension of the `odm` object
#' @export
#' @examples
#' library(sceptredata)
#' directories_to_load <- paste0(
#'  system.file("extdata", package = "sceptredata"),
#'  "/highmoi_example/gem_group_", c(1, 2)
#' )
#' directory_to_write <- tempdir()
#' out_list <- create_odm_from_cellranger(
#'   directories_to_load = directories_to_load,
#'   directory_to_write = directory_to_write,
#' )
#' gene_odm <- out_list$gene
#' # return the dimension, number of rows, and number of columns
#' dim(gene_odm)
#' nrow(gene_odm)
#' ncol(gene_odm)
setMethod(f = "dim", signature = "odm", definition = function(x) x@dimension)
