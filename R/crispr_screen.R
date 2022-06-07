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
