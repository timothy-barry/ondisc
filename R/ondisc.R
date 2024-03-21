utils::globalVariables(c("feature_idx", "vector_idx", "j", "x", "vector_id", "grna_id"))
#' @useDynLib ondisc, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom data.table setkey
#' @examples
#' # initialize odm objects from Cell Ranger output; also, compute the cellwise covariates
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
#'
#' # extract the odm corresponding to the gene modality
#' gene_odm <- out_list$gene
#' gene_odm
#'
#' # obtain dimension information
#' dim(gene_odm)
#' nrow(gene_odm)
#' ncol(gene_odm)
#'
#' # obtain rownames (i.e., the feature IDs)
#' rownames(gene_odm) |> head()
#'
#' # extract row into memory, first by integer and then by string
#' expression_vector_1 <- gene_odm[10,]
#' expression_vector_2 <- gene_odm["ENSG00000135046",]
#'
#' # delete the gene_odm object
#' rm(gene_odm)
#'
#' # reinitialize the gene_odm object
#' gene_odm <- initialize_odm_from_backing_file(
#'   paste0(tempdir(), "/gene.odm")
#' )
#' gene_odm
#' @import Rhdf5lib
"_PACKAGE"
