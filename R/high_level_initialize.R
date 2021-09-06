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
#' @param odm_fp location to write the ondisc matrix to disk
#' @param metadata_fp location to write the me metadata .RDS file. By default, a file called "metadata.rds" stored in the same directory as the backing .odm file.
#' @param barcode_suffixes a vector of suffix that appended to the barcodes to for each input barcode_fp. If NULL, append file index.
#' @param progress (optional; default FALSE) print progress messages?
#' @return A `covariate_ondisc_matrix`.
#' @export
#' @examples
#' # First example: initialize a `covariate_ondisc_matrix`
#' # using simulated expression data; store output in tempdir()
#' tempfile <- create_new_directory()
#' odm_fp <- paste0(tempfile, "/expression_odm")
#' file_locs <- system.file("extdata",package = "ondisc",
#' c("gene_expression.mtx", "genes.tsv", "cell_barcodes.tsv"))
#' names(file_locs) <- c("expressions", "features", "barcodes")
#' expression_data <- create_ondisc_matrix_from_mtx(mtx_fp = file_locs[["expressions"]],
#' barcodes_fp = file_locs[["barcodes"]],
#' features_fp = file_locs[["features"]],
#' odm_fp = odm_fp)
#'
#' # Second example: initialize a `covariate_ondisc_matrix` using simulated
#' # gRNA perturbation data; store in tempdir()
#' tempfile <- create_new_directory()
#' odm_fp <- paste0(tempfile, "/expression_odm")
#' file_locs <- system.file("extdata", package = "ondisc",
#' c("perturbation.mtx", "guides.tsv", "cell_barcodes.tsv"))
#' names(file_locs) <- c("perturbations", "features", "barcodes")
#' perturbation_data <- create_ondisc_matrix_from_mtx(mtx_fp = file_locs[["perturbations"]],
#' barcodes_fp = file_locs[["barcodes"]],
#' features_fp = file_locs[["features"]],
#' odm_fp = odm_fp)
#'
#'
#' # Third example: initialize from a list of .mtx files
#' tempfile <- create_new_directory()
#' odm_fp <- paste0(tempfile, "/expression_odm")
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
#' odm <- create_ondisc_matrix_from_mtx(mtx_fp, barcodes_fp, features_fp, odm_fp)
create_ondisc_matrix_from_mtx <- function(mtx_fp, barcodes_fp, features_fp, odm_fp, metadata_fp = NULL, n_lines_per_chunk = 3e+08, barcode_suffixes = NULL, progress = TRUE) {
  # bag_of_variables is used to store quantities to compute the feature- and cell-covariates.
  # mtx_metadata stores general information about the .mtx file, and features_metadata stores general information about the .tsv files.

  # set the h5 file path
  odm_fp <- append_file_extension(odm_fp, "odm")
  if (file.exists(odm_fp)) stop(paste0("File ", odm_fp, " already exists. Ending function."))

  # generate random ODM id
  odm_id <- sample(seq(0L, .Machine$integer.max), size = 1)

  # Define "bag_of_variables" environment for storing args
  bag_of_variables <- new.env()

  # Extract .mtx metadata; add mtx fp and n_lines per chunk to list
  mtx_metadata <- get_mtx_metadata(mtx_fp)
  mtx_metadata$mtx_fp <- mtx_fp
  mtx_metadata$n_lines_per_chunk <- n_lines_per_chunk
  # update the bag of variables with n_cells, n_features, and n_cells_in_files
  bag_of_variables[[arguments_enum()$n_cells]] <- mtx_metadata$n_cells
  bag_of_variables[[arguments_enum()$n_features]] <- mtx_metadata$n_features
  bag_of_variables[[arguments_enum()$n_cells_in_files]] <- mtx_metadata$n_cells_in_files

  # Extract features.tsv metadata; as a side-effect, if there are MT genes, put the locations of those genes into the bag_of_vars.
  features_metadata <- get_features_metadata(features_fp, bag_of_variables)

  # Initialize the .h5 file on-disk (side-effect)
  initialize_h5_file_on_disk(odm_fp, mtx_metadata, odm_id)

  # Determine which covariates to compute
  covariates <- map_inputs_to_covariates(mtx_metadata, features_metadata)

  # Run core algorithm
  out <- run_core_algo(odm_fp, mtx_metadata, covariates, bag_of_variables, progress) # returns list of 2 covariate matrices (one for cells and one for features); side-effect is to store all information from user input into h5_fp

  # Obtain the string arrays (e.g., feature IDs, feature names, and possibly cell barcodes)
  string_arrays <- get_string_arrays(barcodes_fp, features_fp, features_metadata, barcode_suffixes)

  # initialize metadata ondisc matrix
  odm <- ondisc_matrix(h5_file = odm_fp,
                       logical_mat = mtx_metadata$is_logical,
                       underlying_dimension = c(mtx_metadata$n_features, mtx_metadata$n_cells),
                       feature_ids = string_arrays$feature_ids,
                       feature_names = string_arrays$feature_names,
                       cell_barcodes = string_arrays$cell_barcodes,
                       odm_id = odm_id)

  # initialize the metadata odm
  if (!is.null(barcode_suffixes)) {
    batch <- rep(x = barcode_suffixes, times = cells_metadata$n_cells_in_files) %>% factor()
    out$cell_covariates$batch <- batch
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


#' Create `ondisc_matrix` from h5
#'
#' Creates an `ondisc_matrix` from a list of .h5 files storing single-cell expression data.
#'
#' @param h5_list a list of .h5 files to convert into `ondisc_matrix` format; the .h5 file should have the same features measured on different sets of cells.
#' @param odm_fp location to write the ondisc matrix to disk
#' @param metadata_fp location to write the me metadata .RDS file. By default, a file called "metadata.rds" stored in the same directory as the backing .odm file.
#' @param barcode_suffixes a vector of suffix that appended to the barcodes to for each input barcode_fp. If NULL, append file index.
#' @param progress progress (optional; default FALSE) print progress messages?
#'
#' @return A covariate_ondisc_matrix.
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
#' tempfile <- create_new_directory()
#' odm_fp <- paste0(tempfile, "/expression_odm")
#' # initialize a directory to store ondisc matrix
#' storage_dir <- create_new_directory()
#' # create the ondisc matrix
#' odm_plus_covariates_list <- create_ondisc_matrix_from_h5_list(h5_list, storage_dir)
#' }
create_ondisc_matrix_from_h5_list <- function(h5_list, odm_fp, metadata_fp = NULL, barcode_suffixes = NULL, progress = TRUE) {
  # bag_of_variables is used to store quantities to compute the feature- and cell-covariates.
  # cells_metadata stores general information about the .h5 file and the cells
  # features_metadata stores general information about the features

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
  features_metadata <- get_features_metadata_from_table(features_df, bag_of_variables)

  # Extract cell metadata at the same time
  cells_metadata <- get_h5_cells_metadata(h5_list)
  cells_metadata$h5_list <- h5_list
  cells_metadata$n_features <- nrow(features_df)
  bag_of_variables[[arguments_enum()$n_cells]] <- cells_metadata$n_cells
  bag_of_variables[[arguments_enum()$n_features]] <- cells_metadata$n_features
  bag_of_variables[[arguments_enum()$n_cells_in_files]] <- cells_metadata$n_cells_in_files

  # Initialize the .h5 file on-disk (side-effect)
  initialize_h5_file_on_disk(odm_fp, cells_metadata, odm_id)

  # Determine which covariates to compute
  covariates <- map_inputs_to_covariates(cells_metadata, features_metadata)

  # Run core algorithm
  out <- run_core_algo(odm_fp, cells_metadata, covariates, bag_of_variables, progress) # returns list of 2 covariate matrices (one for cells and one for features); side-effect is to store all information from user input into h5_fp

  # Obtain the string arrays (e.g., feature IDs, feature names, and possibly cell barcodes)
  string_arrays <- list(cell_barcodes = barcodes, feature_ids = features_df$id , feature_names = features_df$name)

  # initialize metadata ondisc matrix
  odm <- ondisc_matrix(h5_file = odm_fp,
                       logical_mat = cells_metadata$is_logical,
                       underlying_dimension = c(cells_metadata$n_features, cells_metadata$n_cells),
                       feature_ids = string_arrays$feature_ids,
                       feature_names = string_arrays$feature_names,
                       cell_barcodes = string_arrays$cell_barcodes,
                       odm_id = odm_id)

  # initialize the metadata odm
  if (!is.null(barcode_suffixes)) {
    batch <- rep(x = barcode_suffixes, times = cells_metadata$n_cells_in_files) %>% factor()
    out$cell_covariates$batch <- batch
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
