#' Create an `ondisc_matrix` from a single .mtx file or a list of .mtx files.
#'
#' Initializes an `ondisc_matrix` from .mtx file(s), a features.tsv file, and barcode.tsv file(s).
#' The number of .mtx files should be the same of the number of the barcodes.tsv files. The list of .mtx files should share the same features.tsv file.
#' Returns an `ondisc_matrix` along with cell-specific and feature-specific covariate matrices.
#'
#' The function can compute the following cell-specific and feature-specific covariates:
#' - cell-specific: (i) total number of features expressed in cell (n_nonzero_cell), (ii) total UMI count (n_umis_cell), and (iii) percentage of UMIs that map to mitochondrial genes (p_mito_cell).
#' - feature-specific: (i) total number of cells in which feature is expressed (n_nonzero_feature), (ii) mean expression of feature across cells (mean_expression_feature), (iii) coefficient of variation of feature expression across cells (coef_of_variation_feature).
#'
#' The function decides which covariates to compute given the input; in general, the function computes the maximum set of covariates possible.
#'
#' @param mtx_fp file path to .mtx file(s) storing the expression data. The .mtx file can represent either an integer matrix or a logical (i.e., binary) matrix. If the .mtx file contains only two columns (after the initial three-column row of metadata), then the .mtx file is assumed to represent a logical matrix.
#' @param barcodes_fp file path to the .tsv file containing the cell barcodes.
#' @param features_fp file path to the features.tsv file. The first column (required) contains the feature IDs (e.g., ENSG00000186092), and the second column (optional) contains the human-readable feature names (e.g., OR4F5). Subsequent columns are discarded.
#' @param n_lines_per_chunk (optional) number of lines in .mtx file to process per chunk. Defaults to 3e+08.
#' @param odm_fp location to write the ondisc matrix to disk.
#' @param metadata_fp (optional; default NULL) location to write the metadata .RDS file. By default, a file called "metadata.rds" stored in the same directory as the backing .odm file.
#' @param barcode_suffixes (optional; default NULL) a vector of suffix that appended to each barcodes in barcodes_fp. The length should be the same as the length of `barcodes_fp`. If NULL, append nothing for a single .mtx input; append file index for a list of .mtx inputs.
#' @param progress (optional; default TRUE) print progress messages?
#' @param comp_level amount of compression to apply (not yet implemented)
#' @return A `covariate_ondisc_matrix`.
#' @export
#'
#' @examples
#' # # Please install the `ondiscdata` package before running the examples.
#' # install.packages("devtools")
#' # devtools::install_github("Katsevich-Lab/ondiscdata")
#'
#' # First example: initialize a `covariate_ondisc_matrix`
#' # using simulated expression data; store output in tempdir()
#' tempfile <- create_new_directory()
#' odm_fp <- paste0(tempfile, "/expression_odm")
#' file_locs <- system.file("extdata/mtx", package = "ondiscdata",
#' c("matrix.mtx", "features.tsv", "barcodes.tsv"))
#' names(file_locs) <- c("expressions", "features", "barcodes")
#' expression_data <- create_ondisc_matrix_from_mtx(mtx_fp = file_locs[["expressions"]],
#' barcodes_fp = file_locs[["barcodes"]],
#' features_fp = file_locs[["features"]],
#' n_lines_per_chunk = 537000,
#' odm_fp = odm_fp)
#'
#' # Second example: initialize from a list of .mtx files
#' tempfile <- create_new_directory()
#' odm_fp <- paste0(tempfile, "/expression_odm")
#' mtx_fp <- system.file("extdata/mtx_list", package = "ondiscdata",
#' c("mtx_dir1/matrix.mtx", "mtx_dir2/matrix.mtx", "mtx_dir3/matrix.mtx"))
#' barcodes_fp <- system.file("extdata/mtx_list", package = "ondiscdata",
#' c("mtx_dir1/barcodes.tsv", "mtx_dir2/barcodes.tsv", "mtx_dir3/barcodes.tsv"))
#' features_fp <- system.file("extdata", "mtx_list/mtx_dir1/features.tsv", package = "ondiscdata")
#' odm <- create_ondisc_matrix_from_mtx(mtx_fp, barcodes_fp, features_fp, odm_fp)
create_ondisc_matrix_from_mtx <- function(mtx_fp, barcodes_fp, features_fp, odm_fp, metadata_fp = NULL, n_lines_per_chunk = 3e+08, barcode_suffixes = NULL, comp_level = 0L, progress = TRUE) {
  # bag_of_variables is used to store quantities to compute the feature- and cell-covariates,
  # general information about the .mtx file, and general information about the .tsv files.

  # set the h5 file path
  odm_fp <- append_file_extension(odm_fp, "odm")
  if (file.exists(odm_fp)) stop(paste0("File ", odm_fp, " already exists. Ending function."))

  # generate random ODM id
  odm_id <- sample(seq(0L, .Machine$integer.max), size = 1)

  # Define "bag_of_variables" environment for storing args
  bag_of_variables <- new.env()

  # Extract .mtx metadata; add mtx fp and n_lines per chunk to list
  list2env(get_mtx_metadata(mtx_fp), bag_of_variables)
  bag_of_variables[["mtx_fp"]] <- mtx_fp
  bag_of_variables[["n_lines_per_chunk"]] <- n_lines_per_chunk

  # Extract features.tsv metadata; as a side-effect, if there are MT genes, put the locations of those genes into the bag_of_vars.
  list2env(get_features_metadata(features_fp, bag_of_variables), bag_of_variables)

  # Initialize the .h5 file on-disk (side-effect)
  initialize_h5_file_on_disk(odm_fp, bag_of_variables, odm_id, comp_level)

  # Determine which covariates to compute
  covariates <- map_inputs_to_covariates(bag_of_variables)

  # Run core algorithm
  out <- run_core_algo(odm_fp, covariates, bag_of_variables, progress) # returns list of 2 covariate matrices (one for cells and one for features); side-effect is to store all information from user input into h5_fp

  # Obtain the string arrays (e.g., feature IDs, feature names, and possibly cell barcodes)
  if (is.null(barcode_suffixes) && length(barcodes_fp) > 1L) {
    barcode_suffixes <- seq(1,length(barcodes_fp))
  }
  string_arrays <- get_string_arrays(barcodes_fp, features_fp, bag_of_variables, barcode_suffixes)


  # initialize metadata ondisc matrix
  odm <- ondisc_matrix(h5_file = odm_fp,
                       logical_mat = bag_of_variables$is_logical,
                       underlying_dimension = c(bag_of_variables$n_features, bag_of_variables$n_cells),
                       feature_ids = string_arrays$feature_ids,
                       feature_names = string_arrays$feature_names,
                       cell_barcodes = string_arrays$cell_barcodes,
                       odm_id = odm_id)

  # initialize the metadata odm
  if (!is.null(barcode_suffixes)) {
    batch <- rep(x = barcode_suffixes, times = bag_of_variables$n_cells_in_files) %>% factor()
    out$cell_covariates$barcode_suffix <- batch
  }
  metadata_odm <- covariate_ondisc_matrix(ondisc_matrix = odm,
                                         cell_covariates = out$cell_covariates,
                                         feature_covariates = out$feature_covariates)

  # save the metadata
  if (is.null(metadata_fp)) metadata_fp <- paste0(dirname(odm_fp), "/metadata.rds")
  save_odm(metadata_odm, metadata_fp)

  # return initialized object
  return(metadata_odm)
}


#' Create `ondisc_matrix` from a single .h5 file or a list of .h5 files.
#'
#' Creates an `ondisc_matrix` from a single .h5 file or a list of .h5 files storing single-cell expression data.
#' The .h5 files should fulfill the following format requirements:
#'   1. There is one and only one dataset named "barcodes";
#'   2. There is one and only one dataset named "data";
#'   3. There is one and only one dataset named "gene_names";
#'   4. There is one and only one dataset named "genes";
#'   5. There is one and only one dataset named "indices";
#'   6. There is one and only one dataset named "indptr";
#'   7. There is one and only one dataset named "shape";
#'   8. They can be in any h5 group;
#'   9. Each .h5 file can fit in memory.
#'
#'
#' The function can compute the following cell-specific and feature-specific covariates:
#' - cell-specific: (i) total number of features expressed in cell (n_nonzero_cell), (ii) total UMI count (n_umis_cell), and (iii) percentage of UMIs that map to mitochondrial genes (p_mito_cell).
#' - feature-specific: (i) total number of cells in which feature is expressed (n_nonzero_feature), (ii) mean expression of feature across cells (mean_expression_feature), (iii) coefficient of variation of feature expression across cells (coef_of_variation_feature).
#'
#' The function decides which covariates to compute given the input; in general, the function computes the maximum set of covariates possible.
#'
#' @param h5_list a single .h5 file or a list of .h5 files to convert into `ondisc_matrix` format; the .h5 file should have the same features measured on different sets of cells.
#' @param odm_fp location to write the ondisc matrix to disk.
#' @param metadata_fp (optional; default NULL) location to write the metadata .RDS file. By default, a file called "metadata.rds" stored in the same directory as the backing .odm file.
#' @param barcode_suffixes (optional; default NULL) a vector of suffix that appended to each  barcodes in input .h5 file(s). The length should be the same as the length of `h5_list`. If NULL, append nothing for a single .h5 input; append file index for a list of .h5 inputs.
#' @param progress progress (optional; default TRUE) print progress messages?
#' @param comp_level amount of compression to apply (not yet implemented)
#'
#' @return A `covariate_ondisc_matrix`.
#' @export
#'
#' @examples
#' #' # Please install the `ondiscdata` package before running the examples.
#' # install.packages("devtools")
#' # devtools::install_github("Katsevich-Lab/ondiscdata")
#'
#' tempfile <- create_new_directory()
#' odm_fp <- paste0(tempfile, "/expression_odm")
#' # get file paths to three .h5 files; these files contain different cells measured on the same genes
#' h5_list <- system.file("extdata/h5_list", package = "ondiscdata",
#' c("batch-1_1.h5", "batch-1_2.h5", "batch_2-1.h5"))
#' # create the ondisc matrix (commented out because of long (~2 min) running time)
#' # odm_plus_covariates_list <- create_ondisc_matrix_from_h5_list(h5_list, odm_fp)
create_ondisc_matrix_from_h5_list <- function(h5_list, odm_fp, metadata_fp = NULL, barcode_suffixes = NULL, comp_level = 0L, progress = TRUE) {
  # bag_of_variables is used to store quantities to compute the feature- and cell-covariates,
  # general information about the .h5 file, the cells and the features

  # set the h5 file path
  odm_fp <- append_file_extension(odm_fp, "odm")
  if (file.exists(odm_fp)) stop(paste0("File ", odm_fp, " already exists. Ending function."))

  # generate random ODM id
  odm_id <- sample(seq(0L, .Machine$integer.max), size = 1)

  # Define "bag_of_variables" environment for storing args
  bag_of_variables <- new.env()

  # Extract barcodes and features_df, compute feature matadata; as a side-effect, if there are MT genes, put the locations of those genes into the bag_of_vars.
  if (is.null(barcode_suffixes) && length(h5_list) > 1L) barcode_suffixes <- seq(1, length(h5_list))
  barcodes_list <- get_h5_barcodes(h5_list, barcode_suffixes)
  barcodes <- unlist(barcodes_list)
  features_df <- get_h5_features(h5_list)
  list2env(get_features_metadata(features_df, bag_of_variables, FALSE), bag_of_variables)

  # Extract cell metadata at the same time
  list2env(get_h5_cells_metadata(h5_list), bag_of_variables)
  bag_of_variables[["h5_list"]] <- h5_list
  bag_of_variables[["n_features"]] <- nrow(features_df)

  # Initialize the .h5 file on-disk (side-effect)
  initialize_h5_file_on_disk(odm_fp, bag_of_variables, odm_id, comp_level)

  # Determine which covariates to compute
  covariates <- map_inputs_to_covariates(bag_of_variables)

  # Run core algorithm
  out <- run_core_algo(odm_fp, covariates, bag_of_variables, progress) # returns list of 2 covariate matrices (one for cells and one for features); side-effect is to store all information from user input into h5_fp

  # Obtain the string arrays (e.g., feature IDs, feature names, and possibly cell barcodes)
  string_arrays <- list(cell_barcodes = barcodes, feature_ids = features_df$id , feature_names = features_df$name)

  # initialize metadata ondisc matrix
  odm <- ondisc_matrix(h5_file = odm_fp,
                       logical_mat = bag_of_variables$is_logical,
                       underlying_dimension = c(bag_of_variables$n_features, bag_of_variables$n_cells),
                       feature_ids = string_arrays$feature_ids,
                       feature_names = string_arrays$feature_names,
                       cell_barcodes = string_arrays$cell_barcodes,
                       odm_id = odm_id)

  # initialize the metadata odm
  if (!is.null(barcode_suffixes)) {
    batch <- rep(x = barcode_suffixes, times = bag_of_variables$n_cells_in_files) %>% factor()
    out$cell_covariates$barcode_suffix <- batch
  }
  metadata_odm <- covariate_ondisc_matrix(ondisc_matrix = odm,
                                          cell_covariates = out$cell_covariates,
                                          feature_covariates = out$feature_covariates)

  # save the metadata
  if (is.null(metadata_fp)) metadata_fp <- paste0(dirname(odm_fp), "/metadata.rds")
  save_odm(metadata_odm, metadata_fp)

  # return initialized object
  return(metadata_odm)
}
