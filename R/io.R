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
#' @return NULL for `save_odm`; a loaded `ondisc` object for `read_odm`
#' @export
#' @examples
#' # Install the `ondiscdata` package to run the examples.
#' # devtools::install_github("timothy-barry/ondiscdata")
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
#' # finally, load the subsetted odm
#' odm_subset <- read_odm(odm_fp, metadata_subset_fp)
save_odm <- function(odm, metadata_fp) {
  if (is(odm, "ondisc_matrix")) {
    odm@h5_file <- NA_character_
  } else if (is(odm, "covariate_ondisc_matrix")) {
    odm@ondisc_matrix@h5_file <- NA_character_
  } else if (is(odm, "multimodal_ondisc_matrix")) {
    stop("Use function `save_multimodal_odm` for multimodal_ondisc_matrix class.")
  } else {
    stop("Object not recognized.")
  }
  metadata_fp <- append_file_extension(metadata_fp, "rds")
  saveRDS(object = odm, file = metadata_fp)
}


#' Read ODM
#' @rdname save_odm
#' @param odm_fp file path to a backing .odm file
#' @export
read_odm <- function(odm_fp, metadata_fp = NULL) {
  if (!file.exists(odm_fp)) stop(paste0("File ", odm_fp, " does not exist."))
  odm_id_backing <- read_integer_vector_hdf5(odm_fp, "odm_id", 1L)

  if (is.null(metadata_fp)) {
    # case 1: load ODM without metadata
    underlying_dimension <- read_integer_vector_hdf5(odm_fp, "dimension", 2L)
    logical_mat <- as.logical(read_integer_vector_hdf5(odm_fp, "logical_mat", 1L))
    odm <- ondisc_matrix(h5_file = odm_fp,
                         logical_mat = logical_mat,
                         underlying_dimension = underlying_dimension,
                         odm_id = odm_id_backing)
  } else {
    # case 2: load ODM with metadata
    if (!file.exists(metadata_fp)) stop(paste0("File ", metadata_fp, " does not exist."))
    odm <- readRDS(metadata_fp)
    if (is(odm, "ondisc_matrix")) {
      odm@h5_file <- odm_fp
      metadata_odm_id <- odm@odm_id
    } else if (is(odm, "covariate_ondisc_matrix")) {
      odm@ondisc_matrix@h5_file <- odm_fp
      metadata_odm_id <- odm@ondisc_matrix@odm_id
    } else {
      stop("Object must be of class ondisc_matrix or covariate_ondisc_matrix.")
    }
    # check that the odm IDs coincide
    if (metadata_odm_id != odm_id_backing) stop("ODM-ID of backing .odm file does not match that of metadata file.")
  }
  return(odm)
}


#' Save or read a `multimodal_ondisc_matrix`
#'
#' Functions to save and read `multimodal_ondisc_matrix` objects.
#'
#' A `multimodal_ondisc_matrix` is associated with two or more backing .odm files. These files are static, read-only files that are never altered or copied.
#' `read_multimodal_odm` constructs a `multimodal_ondisc_matrix` from a given set of backing .odm files and a multimodal metadata RDS file.
#' `save_multimodal_odm` writes a new metadata file to disk (potentially overwriting a previous metadata file of the same name), leaving the backing .odm files unaltered.
#'
#' @param multimodal_odm a `multimodal_ondisc_matrix` object
#' @param multimodal_metadata_fp path to a multimodal metadata RDS file
#'
#' @return NULL for `save_multimodal_odm`; a loaded `multimodal_ondisc_matrix` for `read_multimodal_odm`
#' @export
#'
#' @examples
#' # Install the `ondiscdata` package to run the examples.
#' # devtools::install_github("timothy-barry/ondiscdata")
#'
#' # Initialize a multimodal ODM
#' gene_odm_fp <- system.file("extdata", "odm/gene/matrix.odm", package = "ondiscdata")
#' odm_gene <- read_odm(gene_odm_fp,
#' system.file("extdata", "odm/gene/metadata.rds", package = "ondiscdata"))
#' grna_odm_fp <- system.file("extdata",
#' "odm/grna_assignment/matrix.odm", package = "ondiscdata")
#' odm_grna <- read_odm(grna_odm_fp,
#' system.file("extdata", "odm/grna_assignment/metadata.rds", package = "ondiscdata"))
#' multimodal_odm <- multimodal_ondisc_matrix(list(gene = odm_gene, grna = odm_grna))
#'
#' # subset the multimodal ODM by cell
#' multimodal_odm_sub <- multimodal_odm[,1:10]
#'
#' # save the subsetted multimodal ODM
#' multimodal_metadata_fp <- paste0(create_new_directory(), "/multimodal_metadata.rds")
#' save_multimodal_odm(multimodal_odm_sub, multimodal_metadata_fp)
#'
#' # delete and then read the subsetted multimodal odm
#' rm(multimodal_odm_sub)
#' multimodal_odm_sub <- read_multimodal_odm(c(gene_odm_fp, grna_odm_fp), multimodal_metadata_fp)
save_multimodal_odm <- function(multimodal_odm, multimodal_metadata_fp) {
  modality_names <- names(multimodal_odm@modalities)
  # set h5 fps to NA
  for (modality_name in modality_names) {
    multimodal_odm@modalities[[modality_name]]@ondisc_matrix@h5_file <- NA_character_
  }
  # ensure that there are no ODM ID collisions
  odm_ids <- sapply(multimodal_odm@modalities, function(modality) {
    modality@ondisc_matrix@odm_id
  })
  if (length(unique(odm_ids)) != length(odm_ids)) {
    stop("ODM IDs of different modalities coincide.")
  }
  multimodal_metadata_fp <- append_file_extension(multimodal_metadata_fp, "rds")
  saveRDS(object = multimodal_odm, file = multimodal_metadata_fp)
}


#' Read multimodal ODM
#' @rdname save_multimodal_odm
#' @param odm_fps vector of file paths to backing .odm files
#' @export
read_multimodal_odm <- function(odm_fps, multimodal_metadata_fp) {
  mm_odm <- readRDS(multimodal_metadata_fp)

  # determine which h5 path goes with which modality
  backing_ids <- sapply(X = odm_fps, FUN = function(odm_fp) read_integer_vector_hdf5(odm_fp, "odm_id", 1L))
  memory_ids <- sapply(X = mm_odm@modalities, FUN = function(modality) modality@ondisc_matrix@odm_id)
  # compute the ID intersection
  intersected_ids <- intersect(backing_ids, memory_ids)
  if (length(intersected_ids) == 0) stop("Backing .odm files do not match multimodal_metadata.rds file.")
  odm_fps_ordered <- names(backing_ids)[match(x = intersected_ids, table = backing_ids)]
  modality_names_ordered <- names(memory_ids)[match(x = intersected_ids, table = memory_ids)]
  # assign h5 fp to each modality
  for (i in seq(1, length(odm_fps_ordered))) {
    mm_odm@modalities[[modality_names_ordered[i]]]@ondisc_matrix@h5_file <- odm_fps_ordered[i]
  }
  return(mm_odm)
}
