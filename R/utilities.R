#' setNames
#'
#' My personal setNames function
#'
#' @param object the object to name
#' @param nm the names
#'
#' @return the named object
setNames <- function(object, nm) {
  names(object) <- nm
  object
}
