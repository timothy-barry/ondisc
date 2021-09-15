#' Save or read an `ondisc_matrix`
#'
#' Functions to save and read `ondisc_matrix` objects.
#'
#' The backing .odm file is a static, read-only file that is never altered or copied.
#' `read_odm` constructs an `ondisc_matrix` from a given backing .odm file and (optionally) a metadata file.
#' `save_odm` writes a new metadata file to disk (potentially overwriting a previous metadata file of the same name), leaving the backing .odm file unaltered.
#'
#' @param odm an ondisc_matrix object
#' @param metadata_fp path to metadata RDS file
#'
#' @return NULL
#' @export
#' @examples
#' # Please install the `ondiscdata` package before running the examples.
#' # install.packages("devtools")
#' # devtools::install_github("Katsevich-Lab/ondiscdata")
#'
#' # Load odm from package
#' odm_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata")
#' odm <- read_odm(odm_fp, metadata_fp)
#' # assign a new ondisc_matrix by subsetting the original
#' odm_subset <- odm[1:10,]
#' # save the subsetted odm, and remove it from memory
#' tempfile <- create_new_directory()
#' metadata_subset_fp <- paste0(tempfile, "/metadata_subset.rds")
#' save_odm(odm_subset, metadata_subset_fp)
save_odm <- function(odm, metadata_fp) {
  if (is(odm, "multimodal_ondisc_matrix")) stop("Save function not yet implemented for multimodal_ondisc_matrix class.")

  metadata <- list()
  # if is covariate_ondisc_matrix, save cell-specific and feature-specific covariates, and update odm to the ondisc matrix contained therein.
  if (is(odm, "covariate_ondisc_matrix")) {
    metadata$cell_covariates <- odm@cell_covariates
    metadata$feature_covariates <- odm@feature_covariates
    odm <- odm@ondisc_matrix
  }

  # save aux info
  fields <- c("cell_subset", "feature_subset", "feature_ids", "feature_names", "cell_barcodes", "odm_id")
  for (field in fields) metadata[[field]] <- slot(odm, field)

  # save RDS
  metadata_fp <- append_file_extension(metadata_fp, "rds")
  saveRDS(object = metadata, file = metadata_fp)
}


#' Read ODM
#' @rdname save_odm
#' @param odm_fp file path to a backing .odm file
#' @param metadata_fp metadata.rds file
#' @export
read_odm <- function(odm_fp, metadata_fp = NULL) {
  if (!file.exists(odm_fp)) stop(paste0("File ", odm_fp, " does not exist."))
  # first, obtain underlying dimension and logical_marix boolean
  underlying_dimension <- read_integer_vector_hdf5(odm_fp, "dimension", 2L)
  logical_mat <- as.logical(read_integer_vector_hdf5(odm_fp, "logical_mat", 1L))
  odm_id_backing <- read_integer_vector_hdf5(odm_fp, "odm_id", 1L)
  odm <- ondisc_matrix(h5_file = odm_fp,
                       logical_mat = logical_mat,
                       underlying_dimension = underlying_dimension,
                       odm_id = odm_id_backing)
  # if name is non-null, modify the output futher
  if (!is.null(metadata_fp)) {
    metadata <- readRDS(metadata_fp)
    # verify metadata_odm_id matches backing odm_id
    if (metadata$odm_id != odm_id_backing) {
      stop("ODM-ID of backing .odm file does not match that of metadata file. Try loading a different metadata file.")
    }
    fields <- c("cell_subset", "feature_subset", "feature_ids", "feature_names", "cell_barcodes")
    for (field in fields) slot(odm, field) <- metadata[[field]]
    # if cell covariates and feature covariates are available, initialize the covariate_ondisc_matrix
    if ("cell_covariates" %in% names(metadata) && "feature_covariates" %in% names(metadata)) {
      odm <- covariate_ondisc_matrix(ondisc_matrix = odm,
                                    cell_covariates = metadata$cell_covariates,
                                    feature_covariates = metadata$feature_covariates)
    }
  }
  return(odm)
}
