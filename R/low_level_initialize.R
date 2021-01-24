#' Initialize h5 file on-disk
#'
#' Initialize the on-disk portion on an ondisc_matrix.
#'
#' @param h5_fp file path to the .h5 file to be initialized
#' @param mtx_metadata metadata of the .mtx file
#' @param features_metadata metadata of the features.tsv file
#' @param barcodes_fp file path to the barcodes.tsv file
#' @param features_fp file path to the features.tsv file
#'
#' @return NULL
initialize_h5_file_on_disk <- function(h5_fp, mtx_metadata, features_metadata, barcodes_fp, features_fp) {
  # Create the .h5 file
  rhdf5::h5createFile(h5_fp) %>% invisible()
  # Write metadata
  cell_barcodes <- dplyr::pull(readr::read_tsv(file = barcodes_fp, col_names = FALSE, col_types = "c"))
  rhdf5::h5write(cell_barcodes, h5_fp, "cell_barcodes")
  feature_ids <- read_given_column_of_tsv(col_idx = 1, n_cols = features_metadata$n_cols, tsv_file = features_fp)
  rhdf5::h5write(feature_ids, h5_fp, "feature_ids")
  if (features_metadata$feature_names) {
    feature_names <- read_given_column_of_tsv(col_idx = 2, n_cols = features_metadata$n_cols, tsv_file = features_fp)
    rhdf5::h5write(feature_names, h5_fp, "feature_names")
  }
  rhdf5::h5write(c(mtx_metadata$n_features, mtx_metadata$n_cells), h5_fp, "dimension")
  # Initialize CSC
  rhdf5::h5createDataset(file = h5_fp, dataset = "col_ptr", dims = mtx_metadata$n_cells + 1, storage.mode = "integer", level = 0, chunk = 10) %>% invisible()
  rhdf5::h5createDataset(file = h5_fp, dataset = "row_idxs", dims = mtx_metadata$n_data_points, storage.mode = "integer", level = 0, chunk = 100) %>% invisible()
  if (!mtx_metadata$is_logical) {
  rhdf5::h5createDataset(file = h5_fp, dataset = "data_csc", dims = mtx_metadata$n_data_points, storage.mode = "integer", level = 0, chunk = 100) %>% invisible()
  }
  # Initialize CSR
  rhdf5::h5createDataset(file = h5_fp, dataset = "row_ptr", dims = mtx_metadata$n_features + 1, storage.mode = "integer", level = 0, chunk = 10) %>% invisible()
  rhdf5::h5createDataset(file = h5_fp, dataset = "col_idxs", dims = mtx_metadata$n_data_points, storage.mode = "integer", level = 0, chunk = 100) %>% invisible()
  if (!mtx_metadata$is_logical) {
    rhdf5::h5createDataset(file = h5_fp, dataset = "data_csr", dims = mtx_metadata$n_data_points, storage.mode = "integer", level = 0, chunk = 100) %>% invisible()
  }
  return(invisible())
}
