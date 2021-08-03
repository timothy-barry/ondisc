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
#' # First example: initialize a metadata_ondisc_matrix
#' # using simulated expression data; store output in tempdir()
#' tempfile <- create_new_directory()
#' file_locs <- system.file("extdata",package = "ondisc",
#' c("gene_expression.mtx", "genes.tsv", "cell_barcodes.tsv"))
#' names(file_locs) <- c("expressions", "features", "barcodes")
#' expression_data <- create_ondisc_matrix_from_mtx(mtx_fp = file_locs[["expressions"]],
#' barcodes_fp = file_locs[["barcodes"]],
#' features_fp = file_locs[["features"]],
#' on_disk_dir = tempfile,
#' file_name = "expressions",
#' return_metadata_ondisc_matrix = TRUE)
#'
#' # Second example: initialize a metadata_ondisc_matrix using simulated
#' # gRNA perturbation data; store in tempdir()
#' tempfile <- create_new_directory()
#' file_locs <- system.file("extdata", package = "ondisc",
#' c("perturbation.mtx", "guides.tsv", "cell_barcodes.tsv"))
#' names(file_locs) <- c("perturbations", "features", "barcodes")
#' perturbation_data <- create_ondisc_matrix_from_mtx(mtx_fp = file_locs[["perturbations"]],
#' barcodes_fp = file_locs[["barcodes"]],
#' features_fp = file_locs[["features"]],
#' on_disk_dir = tempfile,
#' file_name = "perturbations",
#' return_metadata_ondisc_matrix = TRUE)
#'
#'
#' # Third example: initialize from a list of .mtx files
#' n_mat <- 5
#' n_row_multi <- 300
#' n_col_multi <- sample(x = seq(100, 300), size = n_mat, replace = TRUE)
#' col_multi_cumsum <- c(0,cumsum(n_col_multi))
#' logical_mat_multi <- vector(mode = "logical", length = n_mat)
#' #' generate the matrices using create_synthetic_data
#' r_mats_plus_data_multi <- vector(mode = "list", length = n_mat)
#' set.seed(1)
#' for (i in seq(1,n_mat)) {
#'   r_mats_plus_data_multi[[i]] <- create_synthetic_data(n_row = n_row_multi,
#'                                                        n_col = n_col_multi[i],
#'                                                        logical_mat = logical_mat_multi[i])
#'   r_mats_plus_data_multi[[i]]$features_df <- r_mats_plus_data_multi[[1]]$features_df
#'   r_mats_plus_data_multi[[i]]$features_fp <- r_mats_plus_data_multi[[1]]$features_fp
#' }
#' mtx_fp <- sapply(X = r_mats_plus_data_multi, function(i) i$matrix_fp)
#' barcodes_fp <- sapply(X = r_mats_plus_data_multi, function(i) i$barcodes_fp)
#' features_fp <- r_mats_plus_data_multi[[1]]$features_fp
#' odm <- create_ondisc_matrix_from_mtx(mtx_fp, barcodes_fp, features_fp,
#' return_metadata_ondisc_matrix = TRUE)
create_ondisc_matrix_from_mtx <- function(mtx_fp, barcodes_fp, features_fp, n_lines_per_chunk = 3e+08, on_disk_dir = NULL, file_name = NULL, return_metadata_ondisc_matrix = FALSE, progress = TRUE) {
  # Define "bag_of_variables" environment for storing args
  bag_of_variables <- new.env()

  # Extract .mtx metadata
  mtx_metadata <- get_mtx_metadata(mtx_fp)
  bag_of_variables[[arguments_enum()$n_cells]] <- mtx_metadata$n_cells
  bag_of_variables[[arguments_enum()$n_features]] <- mtx_metadata$n_features
  bag_of_variables[[arguments_enum()$n_cells_in_files]] <- mtx_metadata$n_cells_in_files

  # Extract features.tsv metadata; as a side-effect, if there are MT genes, put the locations of those genes into the bag_of_vars.
  features_metadata <- get_features_metadata(features_fp, bag_of_variables)
  # Set the on_disk_dir, if necessary
  if (is.null(on_disk_dir)) on_disk_dir <- dirname(mtx_fp[1])

  # Generate a name for the ondisc_matrix .h5 file, if necessary
  if (is.null(file_name)) {
    file_name <- generate_on_disc_matrix_name(on_disk_dir)
  } else {
    if (!grepl(pattern = "*\\.h5$", x = file_name)) file_name <- paste0(file_name, ".h5")
  }
  h5_fp <- file.path(on_disk_dir, file_name)

  # Initialize the .h5 file on-disk (side-effect)
  initialize_h5_file_on_disk(h5_fp, mtx_metadata, features_metadata, barcodes_fp, features_fp, progress, TRUE)

  # Determine which covariates to compute
  covariates <- map_inputs_to_covariates(mtx_metadata, features_metadata)

  # Save is_logical and n_rows_to_skip in variables
  is_logical <- mtx_metadata$is_logical
  n_rows_to_skip <- mtx_metadata$n_rows_to_skip

  # Run core algorithm
  out <- run_core_algo(h5_fp, mtx_fp, is_logical, covariates, bag_of_variables, n_lines_per_chunk, n_rows_to_skip, progress) # returns list of 2 covariate matrices (one for cells and one for features); side-effect is to store all information from user input into h5_fp
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


#' Create `ondisc_matrix` from h5
#'
#' Creates an `ondisc_matrix` from a list of .h5 files storing single-cell expression data.
#'
#' @param h5_list a list of .h5 files to convert into `ondisc_matrix` format; the .h5 file should have the same features measured on different sets of cells.
#' @param on_disk_dir (optional) directory in which to store the on-disk portion of the ondisc_matrix. Defaults to the directory in which the .mtx file is located.
#' @param file_name (optional) name of the file in which to store the .h5 data on-disk. Defaults to ondisc_matrix_x.h5, where x is a unique integer starting at 1.
#' @param return_metadata_ondisc_matrix (optional; default FALSE) return the output as a metadata_ondisc_matrix (instead of a list)
#' @param progress progress (optional; default FALSE) print progress messages?
#'
#' @return A list containing (i) an ondisc_matrix, (ii) a cell-specific covariate matrix, and (iii) a feature-specific covariate matrix; if the parameter return_metadata_ondisc_matrix set to TRUE, converts the list to a metadata_ondisc_matrix before returning.
#' @export
#'
#' @examples
#' \dontrun{
#' # ensure example .h5 data from "crisprdata" package are available
#' devtools::install_github("timothy-barry/crisprdata")
#' # get file paths to three .h5 files; these files contain different cells measured on the same genes
#' f_names <- c("GSM3722728_K562-dCas9-KRAB_5K-sgRNAs_Batch-4_2_filtered_gene_bc_matrices_h5.h5",
#' "GSM3722727_K562-dCas9-KRAB_5K-sgRNAs_Batch-4_1_filtered_gene_bc_matrices_h5.h5",
#' "GSM3722729_K562-dCas9-KRAB_5K-sgRNAs_Batch-1_1_filtered_gene_bc_matrices_h5.h5")
#' h5_list <- system.file("extdata", f_names, package = "crisprdata")
#' # initialize a directory to store ondisc matrix
#' storage_dir <- create_new_directory()
#' # create the ondisc matrix
#' odm_plus_covariates_list <- create_ondisc_matrix_from_h5(h5_list, storage_dir)
#' }
create_ondisc_matrix_from_h5_list <- function(h5_list, on_disk_dir = NULL, file_name = NULL, return_metadata_ondisc_matrix = FALSE, progress = TRUE) {
  # add code here!
}
