#' Create ondisc matrix from R matrix
#'
#' Initializes an ondisc matrix from an R matrix. Returns an `ondisc_matrix` along with cell-specific and feature-specific covariate matrices (or optionally, a `metadata_ondisc_matrix`).
#'
#' @param r_matrix an R matrix. The matrix can be either integer or logical.
#' @param barcodes a character vector giving the cell barcodes.
#' @param features_df a data frame giving the names of the features. The first column (required) contains the feature IDs (e.g., ENSG00000186092), and the second column (optional) contains the human-readable feature names (e.g., OR4F5). Subsequent columns are discarded.
#' @param on_disk_dir directory in which to store the on-disk portion of the ondisc_matrix.
#' @param file_name (optional) name of the file in which to store the .h5 data on-disk. Defaults to ondisc_matrix_x.h5, where x is a unique integer starting at 1.
#' @param return_metadata_ondisc_matrix (optional, default FALSE) return the output as a metadata_ondisc_matrix? FALSE returns a list containing an `ondisc_matrix`, a cell-specific covariate matrix, and a feature-specific covariate matrix. TRUE returns a `metadata_ondisc_matrix.`
#'
#' @return A list containing (i) an ondisc_matrix, (ii) a cell-specific covariate matrix, and (iii) a feature-specific covariate matrix; if the parameter return_metadata_ondisc_matrix set to TRUE, converts the list to a metadata_ondisc_matrix before returning.
#' @export
#'
#' @examples
#' # Examples go here.
create_ondisc_matrix_from_R_matrix <- function(r_matrix, barcodes, features_df, on_disk_dir, file_name = NULL, return_metadata_ondisc_matrix = FALSE) {
  # to do:
  # 1. Write this function.
  # 2. Add one or two examples to the "examples" section of the documentation above.
  # 3. (For later) Write one or two tests to test this function.
}
