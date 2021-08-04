#' Save ODM
#'
#' Saves an ondisc matrix or metadata ondisc matrix.
#'
#' @param odm a metadata_ondisc_matrix object or ondisc_matrix object
#'
#' @return NULL
#' @export
#'
#' @examples
#' tempfile <- create_new_directory()
#' odm_fp <- paste0(tempfile, "/expression_odm")
#' file_locs <- system.file("extdata",package = "ondisc",
#' c("gene_expression.mtx", "genes.tsv", "cell_barcodes.tsv"))
#' names(file_locs) <- c("expressions", "features", "barcodes")
#' expression_data <- create_ondisc_matrix_from_mtx(mtx_fp = file_locs[["expressions"]],
#' barcodes_fp = file_locs[["barcodes"]],
#' features_fp = file_locs[["features"]],
#' odm_fp = odm_fp,
#' return_metadata_ondisc_matrix = TRUE)
#' save_odm(expression_data)
#' expression_data_loaded <- read_odm(odm_fp)
save_odm <- function(odm) {
  aux_info <- list()
  if (is(odm, "multimodal_ondisc_matrix")) {
    stop("Save function not yet implemented for metadata_ondisc_matrix class.
         Save components of class individually.")
  }
  if (is(odm, "metadata_ondisc_matrix")) {
    aux_info$cell_covariates <- odm@cell_covariates
    aux_info$feature_covariates <- odm@feature_covariates
    odm <- odm@ondisc_matrix
  }
  # save aux info
  aux_info$cell_subset <- odm@cell_subset
  aux_info$feature_subset <- odm@feature_subset

  aux_info_fp <- paste0(odm@odm_fp, "/aux_info.rds")
  saveRDS(object = aux_info, file = aux_info_fp)
}


#' Read ODM
#' @rdname save_odm
#' @param odm_fp file path to an ondisc matrix (or metadata ondisc matrix)
read_odm <- function(odm_fp) {
  h5_file <- paste0(odm_fp, "/data.h5")
  aux_info_file <- paste0(odm_fp, "/aux_info.rds")
  underlying_dimension <- as.integer(rhdf5::h5read(file = h5_file, name = "dimension"))
  logical_mat <- as.logical(rhdf5::h5read(file = h5_file, name = "logical_mat"))
  out <- internal_initialize_ondisc_matrix(odm_fp = odm_fp, logical_mat = logical_mat, underlying_dimension = underlying_dimension)
  # read the aux info .rds file
  aux_info <- readRDS(aux_info_file)
  out@cell_subset <- aux_info$cell_subset
  out@feature_subset <- aux_info$feature_subset
  # check if feature covariates and cell covariates are present; if so, coerce to metadata_odm
  if ("cell_covariates" %in% names(aux_info) && "feature_covariates" %in% names(aux_info)) {
    out <- metadata_ondisc_matrix(ondisc_matrix = out,
                           cell_covariates = aux_info$cell_covariates,
                           feature_covariates = aux_info$feature_covariates)
  }
  return(out)
}
