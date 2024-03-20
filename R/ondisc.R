utils::globalVariables(c("feature_idx", "vector_idx", "j", "x", "vector_id", "grna_id"))
#' @useDynLib ondisc, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom data.table setkey
#' @import Rhdf5lib
"_PACKAGE"
