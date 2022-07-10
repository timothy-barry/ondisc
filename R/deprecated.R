read_odm_old <- function(odm_fp, metadata_fp = NULL) {
  if (!file.exists(odm_fp)) stop(paste0("File ", odm_fp, " does not exist."))
  metadata <- if (!is.null(metadata_fp)) readRDS(metadata_fp) else NULL
  read_odm_given_metadata_obj(odm_fp, metadata)
}


save_odm_old <- function(odm, metadata_fp) {
  if (is(odm, "multimodal_ondisc_matrix")) stop("Use function `save_multimodal_odm` for multimodal_ondisc_matrix class.")
  # get list representation
  metadata <- convert_odm_metadata_to_list(odm)
  # save RDS
  metadata_fp <- append_file_extension(metadata_fp, "rds")
  saveRDS(object = metadata, file = metadata_fp)
}


# internal helper function: extracts all relevant odm metadata and puts it into a list
convert_odm_metadata_to_list <- function(odm) {
  metadata <- list()
  # if is covariate_ondisc_matrix, save cell-specific and feature-specific covariates, as well as misc list and post load function info, and update odm to the ondisc matrix contained therein.
  if (is(odm, "covariate_ondisc_matrix")) {
    metadata[["is_covariate_odm"]] <- TRUE
    fields_cov <- c("cell_covariates", "feature_covariates", "misc", "post_load_function", "post_load_function_present")
    for (field in fields_cov) metadata[[field]] <- slot(odm, field)
    odm <- odm@ondisc_matrix
  } else {
    metadata[["is_covariate_odm"]] <- FALSE
  }
  # save aux info
  fields <- c("cell_subset", "feature_subset", "feature_ids", "feature_names", "cell_barcodes", "odm_id")
  for (field in fields) metadata[[field]] <- slot(odm, field)
  return(metadata)
}


# a helper function to load an ODM given a metadata object
read_odm_given_metadata_obj <- function(odm_fp, metadata) {
  underlying_dimension <- read_integer_vector_hdf5(odm_fp, "dimension", 2L)
  logical_mat <- as.logical(read_integer_vector_hdf5(odm_fp, "logical_mat", 1L))
  odm_id_backing <- read_integer_vector_hdf5(odm_fp, "odm_id", 1L)
  odm <- ondisc_matrix(h5_file = odm_fp,
                       logical_mat = logical_mat,
                       underlying_dimension = underlying_dimension,
                       odm_id = odm_id_backing)
  # if name is non-null, modify the output futher
  if (!is.null(metadata)) {
    # verify metadata_odm_id matches backing odm_id
    if (metadata$odm_id != odm_id_backing) {
      stop("ODM-ID of backing .odm file does not match that of metadata file. Try loading a different metadata file.")
    }
    fields <- c("cell_subset", "feature_subset", "feature_ids", "feature_names", "cell_barcodes")
    for (field in fields) slot(odm, field) <- metadata[[field]]
    # if covariate_ondisc_matrix, initialize that
    if (metadata[["is_covariate_odm"]]) {
      odm <- covariate_ondisc_matrix(ondisc_matrix = odm, cell_covariates = metadata[["cell_covariates"]], feature_covariates = metadata[["feature_covariates"]])
      fields_cov <- c("misc", "post_load_function", "post_load_function_present")
      for (field in fields_cov) slot(odm, field) <- metadata[[field]]
    }
  }
  return(odm)
}
