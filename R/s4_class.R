# odm class
#' @export
setClass("odm",
         slots = list(h5_file = "character",
                      dimension = "integer",
                      feature_ids = "character",
                      ptr = "externalptr",
                      integer_id = "integer"))


#' Initialize an ondisc matrix
#'
#' Initializes an object of class `odm` from a backing .odm file.
#'
#' @param odm_file file path to the backing `.odm` file
#'
#' @note The file path `odm_file` should be the fully qualified file path to the backing .odm file. In particular, `odm_file` should not contain the home directory symbol (`~`) or current directory symbol (`.`).
#'
#' @return an object of class `odm`
#' @export
initialize_odm_from_backing_file <- function(odm_file) {
  # 0. checks on odm_file
  odm_file <- expand_tilde(odm_file)
  if (!file.exists(odm_file)) {
    stop(paste0(odm_file, " does not exist."))
  }
  if (!(methods::is(odm_file, "character") && length(odm_file) == 1)) {
    stop(paste0(odm_file, " must be a single string."))
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


#' Load a row of an `odm` into memory
#'
#' The operator `[` can be called on an `odm` to load a row of an `odm` into memory.
#'
#' @param x an object of class `odm`
#' @param i the index of the row to load into memory. `i` can be either an integer (specifying the integer-based index of the row to load into memory) or a string (specifying the row to load into memory by feature ID).
#' @param j not used
#' @param drop not used
#' @export
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


#' dimnames
#'
#' Returns the dimnames of an `odm`.
#'
#' @param x an object of class `odm`
#'
#' @return the dimnames of the `odm`
#' @export
setMethod(f = "dimnames", signature = "odm", definition = function(x) list(x@feature_ids, NULL))


#' dim
#'
#' Returns the dimension (i.e., number of features by number of cells) of an `odm`.
#'
#' @param x an object of class `odm`
#'
#' @return the dimension of the `odm`
#' @export
setMethod(f = "dim", signature = "odm", definition = function(x) x@dimension)
