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
#' # Please install the `ondiscdata` package before running the examples.
#' # install.packages("devtools")
#' # devtools::install_github("Katsevich-Lab/ondiscdata")
#'
#' # EXAMPLE 1: an integer-valued matrix of grna expressions
#' odm_fp <- system.file("extdata", "odm/grna_expression/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata", "odm/grna_expression/metadata.rds", package = "ondiscdata")
#' odm <- read_odm(odm_fp, metadata_fp)
#' indics <- load_thresholded_and_grouped_grna(odm, "grna_group_5")
#' indics <- load_thresholded_and_grouped_grna(odm, c("grna_group_1", "grna_group_11"))
#'
#' # EXAMPLE 2: a logical matrix of grna indicators
#' odm_fp <- system.file("extdata", "odm/grna/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata", "odm/grna/metadata.rds", package = "ondiscdata")
#' odm <- read_odm(odm_fp, metadata_fp)
#' # append grna group information -- three grnas/group
#' odm <- odm |>
#' mutate_feature_covariates(grna_group = rep(paste0("grna_group_", seq(1, nrow(odm)/3)), each = 3))
#' indics <- load_thresholded_and_grouped_grna(odm, c("grna_group_21", "grna_group_6"))
#' # alter grna group information -- one grna/group
#' odm <- odm |>
#' mutate_feature_covariates(grna_group = paste0("grna_group_", seq(1, nrow(odm))))
#' indics <- load_thresholded_and_grouped_grna(odm, "grna_group_2")
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
