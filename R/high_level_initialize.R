#' Create an `ondisc_matrix` from a .mtx file.
#'
#' Initializes an `ondisc_matrix` from a .mtx file, a features.tsv file, and a barcodes.tsv file. Returns an `ondisc_matrix` along with cell-specific and feature-specific covariate matrices.
#'
#' The function can compute the following cell-specific and feature-specific covariates:
#' - cell-specific: (i) total number of features expressed in cell (n_nonzero_cell), (ii) total UMI count (n_umis_cell), and (iii) percentage of UMIs that map to mitochondrial genes (p_mito_cell).
#' - feature-specific: (i) total number of cells in which feature is expressed (n_nonzero_feature), (ii) mean expression of feature across cells (mean_expression_feature), (iii) coefficient of variation of feature expression across cells (coef_of_variation_feature).
#'
#' The function decides which covariates to compute given the input; in general, the function computes the maximum set of covariates possible.
#'
#' @param mtx_fp file path to a .mtx file storing the expression data. The .mtx file can represent either an integer matrix or a logical (i.e., binary) matrix. If the .mtx file contains only two columns (after the initial three-column row of metadata), then the .mtx file is assumed to represent a logical matrix.
#' @param barcodes_fp file path to the .tsv file containing the cell barcodes.
#' @param features_fp file path to the features.tsv file. The first column (required) contains the feature IDs (e.g., ENSG00000186092), and the second column (optional) contains the human-readable feature names (e.g., OR4F5). Subsequent columns are discarded.
#' @param n_lines_per_chunk (optional) number of lines in .mtx file to process per chunk. Defaults to 3e+08.
#' @param on_disk_dir (optional) directory in which to store the on-disk portion of the ondisc_matrix. Defaults to the directory in which the .mtx file is located.
#' @param file_name (optional) name of the file in which to store the .h5 data on-disk. Defaults to ondisc_matrix_x.h5, where x is a unique integer starting at 1.
#' @param return_metadata_ondisc_matrix (optional) return the output as a metadata_ondisc_matrix (instead of a list)? Defaults to FALSE.
#' @param progress (optional; default FALSE) print progress messages?
#' @return A list containing (i) an ondisc_matrix, (ii) a cell-specific covariate matrix, and (iii) a feature-specific covariate matrix; if the parameter return_metadata_ondisc_matrix set to TRUE, converts the list to a metadata_ondisc_matrix before returning.
#' @export
#' @examples
#' \dontrun{
#' # First example: initialize a metadata_ondisc_matrix
#' # using simulated expression data; store output in tempdir()
#' file_locs <- system.file("extdata",package = "ondisc",
#' c("gene_expression.mtx", "genes.tsv", "cell_barcodes.tsv"))
#' names(file_locs) <- c("expressions", "features", "barcodes")
#' expression_data <- create_ondisc_matrix_from_mtx(mtx_fp = file_locs[["expressions"]],
#' barcodes_fp = file_locs[["barcodes"]],
#' features_fp = file_locs[["features"]],
#' on_disk_dir = tempdir(),
#' file_name = "expressions",
#' return_metadata_ondisc_matrix = TRUE)
#' saveRDS(object = expression_data, file = paste0(tempdir(), "/expressions.rds"))
#'
#' # Second example: initialize a metadata_ondisc_matrix using simulated
#' # gRNA perturbation data; store in tempdir()
#' file_locs <- system.file("extdata", package = "ondisc",
#' c("perturbation.mtx", "guides.tsv", "cell_barcodes.tsv"))
#' names(file_locs) <- c("perturbations", "features", "barcodes")
#' perturbation_data <- create_ondisc_matrix_from_mtx(mtx_fp = file_locs[["perturbations"]],
#' barcodes_fp = file_locs[["barcodes"]],
#' features_fp = file_locs[["features"]],
#' on_disk_dir = tempdir(),
#' file_name = "perturbations",
#' return_metadata_ondisc_matrix = TRUE)
#' saveRDS(object = perturbation_data, file = paste0(tempdir(), "/perturbations.rds"))
#' }
create_ondisc_matrix_from_mtx <- function(mtx_fp, barcodes_fp, features_fp, n_lines_per_chunk = 3e+08, on_disk_dir = NULL, file_name = NULL, return_metadata_ondisc_matrix = FALSE, progress = TRUE) {
  # Define "bag_of_variables" environment for storing args
  bag_of_variables <- new.env()

  # extract .mtx metadata
  mtx_metadata <- get_mtx_metadata(mtx_fp)
  bag_of_variables[[arguments_enum()$n_cells]] <- mtx_metadata$n_cells
  bag_of_variables[[arguments_enum()$n_features]] <- mtx_metadata$n_features

  # extract features.tsv metadata; as a side-effect, if there are MT genes, put the locations of those genes into the bag_of_vars.
  features_metadata <- get_features_metadata(features_fp, bag_of_variables)
  # set the on_disk_dir, if necessary
  if (is.null(on_disk_dir)) on_disk_dir <- gsub(pattern = '/[^/]*$', replacement = "", x = mtx_fp)

  # Generate a name for the ondisc_matrix .h5 file, if necessary
  if (is.null(file_name)) {
    file_name <- generate_on_disc_matrix_name(on_disk_dir)
  } else {
    if (!grepl(pattern = "*.h5$", x = file_name)) file_name <- paste0(file_name, ".h5")
  }
  h5_fp <- paste0(on_disk_dir, "/", file_name)

  # Initialize the .h5 file on-disk (side-effect)
  initialize_h5_file_on_disk(h5_fp, mtx_metadata, features_metadata, barcodes_fp, features_fp, progress)

  # Determine which covariates to compute
  covariates <- map_inputs_to_covariates(mtx_metadata, features_metadata)

  # Get number of elements to load per chunk
  is_logical <- mtx_metadata$is_logical
  n_rows_to_skip <- mtx_metadata$n_rows_to_skip

  # Run core algorithm
  out <- run_core_mtx_algo(h5_fp, mtx_fp, is_logical, covariates, bag_of_variables, n_lines_per_chunk, n_rows_to_skip, progress)
  odm <- internal_initialize_ondisc_matrix(h5_file = h5_fp, logical_mat = is_logical, underlying_dimension = c(mtx_metadata$n_features, mtx_metadata$n_cells))
  out$ondisc_matrix <- odm
  if (return_metadata_ondisc_matrix) {
    out <- metadata_ondisc_matrix(ondisc_matrix = out$ondisc_matrix,
                                   cell_covariates = out$cell_covariates,
                                   feature_covariates = out$feature_covariates)
  }
  return(out)
}


#' internal initialize ondisc_matrix
#'
#' An internal function for initializing an ondisc_matrix
#'
#' @param h5_file an h5 file
#' @param logical_mat logical value indicating whether the matrix is logical
#' @param underlying_dimension length-two integer vector giving the dimension of the underlying matrix
#'
#' @return a correctly-initialized ondisc_matrix
#' @noRd
internal_initialize_ondisc_matrix <- function(h5_file, logical_mat, underlying_dimension) {
  out <- new(Class = "ondisc_matrix")
  out@h5_file <- h5_file
  out@logical_mat <- logical_mat
  out@underlying_dimension <- underlying_dimension
  return(out)
}
