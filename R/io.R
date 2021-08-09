#' Save ODM
#'
#' Saves an ondisc matrix (or metadata ondisc matrix).
#'
#' This function does not change the underlying .odm expression file. Instead, it creates a new metadata.rds file representing a new "version" of the underlying expression matrix.
#'
#' @param odm a metadata_ondisc_matrix object or ondisc_matrix object
#' @param metadata_fp file path to metadata.rds object
#'
#' @return NULL
#' @export
save_odm <- function(odm, metadata_fp) {
  if (is(odm, "multimodal_ondisc_matrix")) stop("Save function not yet implemented for multimodal_ondisc_matrix class.")

  aux_info <- list()
  # if is metadata_ondisc_matrix, save cell-specific and feature-specific covariates, and update odm to the ondisc matrix contained therein.
  if (is(odm, "metadata_ondisc_matrix")) {
    aux_info$cell_covariates <- odm@cell_covariates
    aux_info$feature_covariates <- odm@feature_covariates
    odm <- odm@ondisc_matrix
  }

  # save aux info
  fields <- c("cell_subset", "feature_subset", "feature_ids", "feature_names", "cell_barcodes")
  for (field in fields) aux_info[[field]] <- slot(odm, field)

  # save RDS
  metadata_fp <- append_file_extension(metadata_fp, "rds")
  saveRDS(object = aux_info, file = metadata_fp)
}


#' Read ODM
#' @rdname save_odm
#' @export
#' @param odm_fp file path to a .odm file
read_odm <- function(odm_fp, metadata_fp = NULL) {
  # first, obtain underlying dimension and logical_marix boolean
  underlying_dimension <- read_integer_vector_hdf5(odm_fp, "dimension", 2L)
  logical_mat <- as.logical(read_integer_vector_hdf5(odm_fp, "logical_mat", 1L))
  odm <- ondisc_matrix(h5_file = odm_fp,
                       logical_mat = logical_mat,
                       underlying_dimension = underlying_dimension)
  # if name is non-null, modify the output futher
  if (!is.null(metadata_fp)) {
    metadata <- readRDS(metadata_fp)
    fields <- c("cell_subset", "feature_subset", "feature_ids", "feature_names", "cell_barcodes")
    for (field in fields) slot(odm, field) <- metadata[[field]]
    # if cell covariates and feature covariates are available, initialize the metadata_ondisc_matrix
    if ("cell_covariates" %in% names(metadata) && "feature_covariates" %in% names(metadata)) {
      odm <- metadata_ondisc_matrix(ondisc_matrix = odm,
                                    cell_covariates = metadata$cell_covariates,
                                    feature_covariates = metadata$feature_covariates)
    }
  }
  return(odm)
}
