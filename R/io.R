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


#' Read ODM
#' @rdname save_odm
#' @param odm_fp file path to a backing .odm file
#' @param metadata_fp metadata.rds file
#' @export
read_odm <- function(odm_fp, metadata_fp = NULL) {
  if (!file.exists(odm_fp)) stop(paste0("File ", odm_fp, " does not exist."))
  metadata <- if (!is.null(metadata_fp)) readRDS(metadata_fp) else NULL
  read_odm_given_metadata_obj(odm_fp, metadata)
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


#' Save a multimodal `ondisc_matrix`
#'
#' @param multimodal_odm a `multimodal_ondisc_matrix` object
#' @param multimodal_metadata_fp path to multimodal metadata RDS file
#'
#' @return NULL
#' @export
#'
#' @examples
#' # Please install the `ondiscdata` package before running the examples.
#' # install.packages("devtools")
#' # devtools::install_github("Katsevich-Lab/ondiscdata")
#'
#' # Initialize a multimodal ODM
#' gene_odm_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' odm_gene <- read_odm(gene_odm_fp, system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata"))
#' gRNA_odm_fp <- system.file("extdata", "odm/gRNA/matrix.odm", package = "ondiscdata")
#' odm_gRNA <- read_odm(gRNA_odm_fp, system.file("extdata", "odm/gRNA/metadata.rds", package = "ondiscdata"))
#' multimodal_odm <- multimodal_ondisc_matrix(list(gene = odm_gene, gRNA = odm_gRNA))
#' multimodal_metadata_fp <- paste0(create_new_directory(), "/multimodal_metadata.rds")
#' save_multimodal_odm(multimodal_odm, multimodal_metadata_fp)
save_multimodal_odm <- function(multimodal_odm, multimodal_metadata_fp) {
  # convert the individual modalities to list form
  modality_lists <- lapply(X = multimodal_odm@modalities, FUN = function(modality) {
    convert_odm_metadata_to_list(modality)
  })
  # create the multimodal metadata list
  metadata <- list(global_cell_covariates = multimodal_odm@global_cell_covariates,
                   modality_lists = modality_lists)
  # save RDS
  multimodal_metadata_fp <- append_file_extension(multimodal_metadata_fp, "rds")
  saveRDS(object = metadata, file = multimodal_metadata_fp)
}


# odm_fps <- c(gene_odm_fp, gRNA_odm_fp)
# multimodal_metadata_fp
read_multimodal_odm <- function(odm_fps, multimodal_metadata_fp) {
  # Load the metadata
  metadata_list <- readRDS(multimodal_metadata_fp)
  # get the global covariate matrix
  global_cell_covariates <- metadata_list$global_cell_covariates
  # determine the order in which to initialize the odms
  backing_ids <- sapply(X = odm_fps, FUN = function(odm_fp) {
    read_integer_vector_hdf5(odm_fp, "odm_id", 1L)
  }) |> sort()
  memory_ids <- sapply(X = metadata$modality_lists, FUN = function(modality) {
    modality$odm_id
  }) |> sort()
  # obtain the ordered odm_fps and modality names
  odm_fps_ordered <- names(backing_ids)
  modality_names_ordered <- names(memory_ids)
  # load the modalities
  modality_list <- lapply(X = seq(1, length(odm_fps_ordered)), FUN = function(i) {
    read_odm_given_metadata_obj(odm_fps_ordered[i],
                                metadata_list$modality_lists[[modality_names_ordered[i]]])
  }) |> setNames(modality_names_ordered)
  # finally, create the multimodal ondisc matrix
  out <- new(Class = "multimodal_ondisc_matrix")
  out@modalities <- covariate_ondisc_matrix_list
  out@global_cell_covariates <- global_cell_covariates
  return(out)
}
