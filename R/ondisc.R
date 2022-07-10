utils::globalVariables(c("gene_idx", "expression", ".", "cells_metadata", "p_mito", "n_umis", "n_nonzero", "p_exp", "coef_of_variation", ".intercept", "in_range", "target", "grna_group"))

#' ondisc: A package for out-of-memory computing on single-cell data
#'
#' Single-cell datasets are large and are growing in size as sequencing costs drop. The ondisc package is designed to facilitate large-scale computing on single-cell expression data by providing access to expression matrices out-of-memory. ondisc is functional (i.e., all objects are persistent) and efficient (i.e., all algorithms are theoretically optimal in time).
#' @useDynLib ondisc, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr %>%
#' @import methods
#' @import Matrix
#' @docType package
#'
#' @name ondisc
NULL
