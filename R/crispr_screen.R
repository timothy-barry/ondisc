#' Convert a list of cell-to-grna assignments to sparse ODM matrix format
#'
#' @param cell_barcodes a vector of cell barcodes
#' @param grna_ids a vector of grna IDs
#' @param grna_assignment_list a list giving the cell-to-grna assignments. The `i`th entry of the list should be a character vector giving the grnas present in cell `i`. If no grnas are present in cell `i`, then an empty string ("") should be supplied.
#' @param odm_fp file path to the backing .odm file to write
#' @param metadata_fp file path to the metadata.rds file to write
#' @param features_metadata_df (optional) a data frame to concatenate to the `feature_covariate_matrix` of the created ODM.
#'
#' @return a logical ODM representing the cell-to-grna assignments. As a side-effect, writes the created ODM to disk
#' @export
#'
#' @examples
#' set.seed(5)
#' n_cells <- 50
#' n_grnas <- 30
#' cell_barcodes <- paste0("cell_", seq(1, n_cells))
#' grna_ids <- paste0("grna_", seq(1, n_grnas))
#' grna_assignment_list <- replicate(n_cells,
#' sample(x = grna_ids, size = sample(x = seq(0, 3), size = 1), replace = FALSE),
#' FALSE) |> lapply(FUN = function(x) if (length(x) == 0) "" else x)
#' features_metadata_df <- data.frame(target_type = sample(x = c("gene", "non-targeting"),
#' size = n_grnas, replace = TRUE),
#' target = paste0("target", seq(1, n_grnas)))
#' odm_fp <- paste0(tempdir(), "/matrix.odm")
#' metadata_fp <- paste0(tempdir(), "/metadata.rds")
#' odm <- convert_assign_list_to_sparse_odm(cell_barcodes, grna_ids, grna_assignment_list,
#' odm_fp, metadata_fp, features_metadata_df)
convert_assign_list_to_sparse_odm <- function(cell_barcodes, grna_ids, grna_assignment_list, odm_fp, metadata_fp, features_metadata_df = NULL) {
  n_cells <- length(cell_barcodes)
  n_grnas <- length(grna_ids)
  # replace "" with character(0) in cell barcodes
  grna_assignment_list <- lapply(grna_assignment_list, function(entry) if (identical(entry, "")) character(0) else entry)
  # get the cell coordinates
  cell_coordinates <- lapply(X = seq(1, n_cells), FUN = function(i) {
    grna_assignment_list[[i]]
    rep(i, length(grna_assignment_list[[i]]))
  }) |> unlist()
  # get the grna coordinates
  grna_coorindates <- sapply(unlist(grna_assignment_list), function(grna_id) which(grna_id == grna_ids))
  # construct the sparse logical matrix
  grna_assignment_m <- as.matrix(Matrix::sparseMatrix(i = grna_coorindates,
                                                      j = cell_coordinates,
                                                      dims = c(n_grnas, n_cells)))
  # write the ODM
  ret <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = grna_assignment_m,
                                                    barcodes = cell_barcodes,
                                                    features_df = data.frame(grna_ids),
                                                    odm_fp = odm_fp,
                                                    metadata_fp = metadata_fp)
  if (!is.null(features_metadata_df)) {
    ret <- ret |> ondisc::mutate_feature_covariates(features_metadata_df)
    ondisc::save_odm(odm = ret, metadata_fp = metadata_fp)
  }
  return(ret)
}


#' Load thresholded and grouped grna
#'
#' Loads data from an ondisc matrix into memory after thresholding and grouping grnas.
#'
#' @param covariate_odm a grna-by-cell covariate ondisc matrix. The matrix can be either an integer-valued expression matrix or a logical matrix of grna-to-cell assignments
#' @param grna_group the grna group (or vector of grna groups) to load
#' @param threshold the threshold to use for assigning grnas to cells (ignored if `covariate_odm` is a logical matrix)
#' @param grna_group_name the name of the column of the feature covariate matrix containing the grna group information
#'
#' @return a matrix of grouped grna-to-cell assignments
#' @export
#'
#' @examples
#' # Install the `ondiscdata` package to run the examples.
#' # devtools::install_github("timothy-barry/ondiscdata")
#'
#' # EXAMPLE 1: an integer-valued matrix of grna expressions
#' odm_fp <- system.file("extdata", "odm/grna_expression/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata",
#' "odm/grna_expression/metadata.rds", package = "ondiscdata")
#' odm <- read_odm(odm_fp, metadata_fp)
#' grps_to_extract <- odm |> get_feature_covariates() |> dplyr::pull(grna_group) |> sample(2)
#' indics <- load_thresholded_and_grouped_grna(odm, grps_to_extract)
#'
#' # EXAMPLE 2: a logical matrix of grna indicators
#' odm_fp <- system.file("extdata", "odm/grna_assignment/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata",
#' "odm/grna_assignment/metadata.rds", package = "ondiscdata")
#' odm <- read_odm(odm_fp, metadata_fp)
#' indics <- load_thresholded_and_grouped_grna(odm, grps_to_extract)
load_thresholded_and_grouped_grna <- function(covariate_odm, grna_group, threshold = 3, grna_group_name = "grna_group") {
  grna_group_vect <- get_feature_covariates(covariate_odm)[[grna_group_name]]
  out <- t(sapply(X = grna_group, FUN = function(curr_grna_group) {
    m <- covariate_odm[[which(curr_grna_group == grna_group_vect),]]
    # first, perform threshold
    m_thresh <- if (!covariate_odm@ondisc_matrix@logical_mat) m >= threshold else m
    # next, if applicable, perform the grouping
    if (nrow(m_thresh) >= 2) Matrix::colSums(m_thresh) >= 1 else as.matrix(m_thresh)
  }))
  return(out)
}


#' Thin multimodal ODM
#'
#' "Thins" a multimodal ODM representing a single-cell CRISPR screen experiment by removing superflous and duplicate information.
#'
#' @param multimodal_odm a multimodal ondisc matrix
#' @param grna_modality_name the name of the grna modality
#' @param gene_modality_name the name of the gene modality
#'
#' @return a "thinned" multimodal ondisc matrix
thin_multimodal_odm <- function(multimodal_odm, grna_modality_name = "grna", response_modality_name = "response", grna_group_column_name = "grna_group") {
  row.names(multimodal_odm@global_cell_covariates) <- NULL

  multimodal_odm@modalities <- multimodal_odm@modalities[c(grna_modality_name, response_modality_name)]
  # clear out (i) cell covariates df, (ii) misc list, (iii) feature names, (iv) cell barcodes
  for (modality in c(grna_modality_name, response_modality_name)) {
    modality_dim <- dim(multimodal_odm@modalities[[modality]])
    multimodal_odm@modalities[[modality]]@cell_covariates <- as.data.frame(matrix(nrow = modality_dim[2], ncol = 0))
    multimodal_odm@modalities[[modality]]@misc <- list()
    multimodal_odm@modalities[[modality]]@ondisc_matrix@feature_names <- ""
    multimodal_odm@modalities[[modality]]@ondisc_matrix@cell_barcodes <- ""
  }

  # clear out feature df of response modality
  multimodal_odm@modalities[[response_modality_name]]@feature_covariates <- data.frame()
  # clear out columns of feature df of grna modality that are not "grna_group"
  multimodal_odm@modalities[[grna_modality_name]]@feature_covariates <-
    multimodal_odm@modalities[[grna_modality_name]]@feature_covariates |> dplyr::select(!!grna_group_column_name)

  return(multimodal_odm)
}
