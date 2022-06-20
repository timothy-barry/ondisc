#' Convert a list of cell-to-gRNA assignments to sparse ODM matrix format
#'
#' @param cell_barcodes a vector of cell barcodes
#' @param gRNA_ids a vector of gRNA IDs
#' @param gRNA_assignment_list a list giving the cell-to-gRNA assignments. The `i`th entry of the list should be a character vector giving the gRNAs present in cell `i`. If no gRNAs are present in cell `i`, then an empty string ("") should be supplied.
#' @param odm_fp file path to the backing .odm file to write
#' @param metadata_fp file path to the metadata.rds file to write
#' @param features_metadata_df (optional) a data frame to concatenate to the `feature_covariate_matrix` of the created ODM.
#'
#' @return a logical ODM representing the cell-to-gRNA assignments. As a side-effect, writes the created ODM to disk
#' @export
#'
#' @examples
#' set.seed(5)
#' n_cells <- 50
#' n_gRNAs <- 30
#' cell_barcodes <- paste0("cell_", seq(1, n_cells))
#' gRNA_ids <- paste0("gRNA_", seq(1, n_gRNAs))
#' gRNA_assignment_list <- replicate(n_cells,
#' sample(x = gRNA_ids, size = sample(x = seq(0, 3), size = 1), replace = FALSE),
#' FALSE) |> lapply(FUN = function(x) if (length(x) == 0) "" else x)
#' features_metadata_df <- data.frame(target_type = sample(x = c("gene", "non-targeting"),
#' size = n_gRNAs, replace = TRUE),
#' target = paste0("target", seq(1, n_gRNAs)))
#' odm_fp <- paste0(tempdir(), "/matrix.odm")
#' metadata_fp <- paste0(tempdir(), "/metadata.rds")
#' odm <- convert_assign_list_to_sparse_odm(cell_barcodes, gRNA_ids, gRNA_assignment_list,
#' odm_fp, metadata_fp, features_metadata_df)
convert_assign_list_to_sparse_odm <- function(cell_barcodes, gRNA_ids, gRNA_assignment_list, odm_fp, metadata_fp, features_metadata_df = NULL) {
  n_cells <- length(cell_barcodes)
  n_gRNAs <- length(gRNA_ids)
  # replace "" with character(0) in cell barcodes
  gRNA_assignment_list <- lapply(gRNA_assignment_list, function(entry) if (identical(entry, "")) character(0) else entry)
  # get the cell coordinates
  cell_coordinates <- lapply(X = seq(1, n_cells), FUN = function(i) {
    gRNA_assignment_list[[i]]
    rep(i, length(gRNA_assignment_list[[i]]))
  }) |> unlist()
  # get the gRNA coordinates
  gRNA_coorindates <- sapply(unlist(gRNA_assignment_list), function(gRNA_id) which(gRNA_id == gRNA_ids))
  # construct the sparse logical matrix
  gRNA_assignment_m <- as.matrix(Matrix::sparseMatrix(i = gRNA_coorindates,
                                                      j = cell_coordinates,
                                                      dims = c(n_gRNAs, n_cells)))
  # write the ODM
  ret <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = gRNA_assignment_m,
                                                    barcodes = cell_barcodes,
                                                    features_df = data.frame(gRNA_ids),
                                                    odm_fp = odm_fp,
                                                    metadata_fp = metadata_fp)
  if (!is.null(features_metadata_df)) {
    ret <- ret |> ondisc::mutate_feature_covariates(features_metadata_df)
    ondisc::save_odm(odm = ret, metadata_fp = metadata_fp)
  }
  return(ret)
}



#' Load thresholded and grouped gRNA
#'
#' Loads data from a covariate ondisc matrix into memory by thresholding and grouping gRNAs.
#'
#' @param covariate_odm a gRNA-by-cell covariate ondisc matrix. The matrix can be either an integer-valued expression matrix or a logical matrix of gRNA-to-cell assignments
#' @param gRNA_group the gRNA group (or vector of gRNA groups) to load
#' @param threshold the threshold to use for assigning gRNAs to cells (ignored if `covariate_odm` is a logical matrix)
#' @param gRNA_group_name the name of the column of the feature covariate matrix containing the gRNA group information
#'
#' @return a matrix of grouped gRNA-to-cell assignments
#' @export
#'
#' @examples
#' # Please install the `ondiscdata` package before running the examples.
#' # install.packages("devtools")
#' # devtools::install_github("Katsevich-Lab/ondiscdata")
#'
#' # EXAMPLE 1: an integer-valued matrix of gRNA expressions
#' odm_fp <- system.file("extdata", "odm/gRNA_expression/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata", "odm/gRNA_expression/metadata.rds", package = "ondiscdata")
#' odm <- read_odm(odm_fp, metadata_fp)
#' indics <- load_thresholded_and_grouped_gRNA(odm, "gRNA_group_5")
#' indics <- load_thresholded_and_grouped_gRNA(odm, c("gRNA_group_1", "gRNA_group_11"))
#'
#' # EXAMPLE 2: a logical matrix of gRNA indicators
#' odm_fp <- system.file("extdata", "odm/gRNA/matrix.odm", package = "ondiscdata")
#' metadata_fp <- system.file("extdata", "odm/gRNA/metadata.rds", package = "ondiscdata")
#' odm <- read_odm(odm_fp, metadata_fp)
#' # append gRNA group information -- three gRNAs/group
#' odm <- odm |>
#' mutate_feature_covariates(gRNA_group = rep(paste0("gRNA_group_", seq(1, nrow(odm)/3)), each = 3))
#' indics <- load_thresholded_and_grouped_gRNA(covariate_odm, c("gRNA_group_21", "gRNA_group_6"))
#' # alter gRNA group information -- one gRNA/group
#' odm <- odm |>
#' mutate_feature_covariates(gRNA_group = paste0("gRNA_group_", seq(1, nrow(odm))))
#' indics <- load_thresholded_and_grouped_gRNA(odm, "gRNA_group_2")
load_thresholded_and_grouped_gRNA <- function(covariate_odm, gRNA_group, threshold = 3, gRNA_group_name = "gRNA_group") {
  gRNA_group_vect <- get_feature_covariates(covariate_odm)[[gRNA_group_name]]
  out <- t(sapply(X = gRNA_group, FUN = function(curr_gRNA_group) {
    m <- covariate_odm[[which(curr_gRNA_group == gRNA_group_vect),]]
    # first, perform threshold
    m_thresh <- if (!covariate_odm@ondisc_matrix@logical_mat) m >= threshold else m
    # next, if applicable, perform the grouping
    if (nrow(m_thresh) >= 2) Matrix::colSums(m_thresh) >= 1 else as.matrix(m_thresh)
  }))
  return(out)
}
