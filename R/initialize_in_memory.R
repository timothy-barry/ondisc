#' Create an `ondisc` matrix from an R matrix
#'
#' Initializes an `ondisc_matrix` from an R matrix. Returns an `ondisc_matrix` alongside cell-specific and feature-specific covariate matrices (or optionally, a `metadata_ondisc_matrix`).
#'
#' This function computes the following cell-specific and feature-specific covariates:
#' - cell-specific: (i) total number of features expressed in cell (n_nonzero_cell), (ii) total UMI count (n_umis_cell), and (iii) percentage of UMIs that map to mitochondrial genes (p_mito_cell).
#' - feature-specific: (i) total number of cells in which feature is expressed (n_nonzero_feature), (ii) mean expression of feature across cells (mean_expression_feature), (iii) coefficient of variation of feature expression across cells (coef_of_variation_feature).
#'
#' @param r_matrix an R matrix. The matrix can be represented as a standard (dense) R matrix or a sparse matrix of class "dgTMatrix". Integer and logical (i.e., boolean) matrices are allowed.
#' @param barcodes a character vector giving the cell barcodes.
#' @param features_df a data frame giving the names of the features. The first column (required) contains the feature IDs (e.g., ENSG00000186092), and the second column (optional) contains the human-readable feature names (e.g., OR4F5). Subsequent columns are discarded. Gene names starting with "MT-" are assumed to be mitochondrial genes and will be used to compute the p_mito covariate.
#' @param on_disk_dir directory in which to store the .h5 file.
#' @param file_name (optional) name of the file in which to store the .h5 data on-disk. Defaults to ondisc_matrix_x.h5, where x is a unique integer starting at 1.
#' @param return_metadata_ondisc_matrix (optional, default FALSE) return the output as a metadata_ondisc_matrix? FALSE returns a list containing an `ondisc_matrix`, a cell-specific covariate matrix, and a feature-specific covariate matrix. TRUE returns a `metadata_ondisc_matrix.`
#'
#' @return A list containing (i) an ondisc_matrix, (ii) a cell-specific covariate matrix, and (iii) a feature-specific covariate matrix; if the parameter return_metadata_ondisc_matrix set to TRUE, converts the list to a `metadata_ondisc_matrix` before returning.
#' @export
#'
#' @examples
#' ##################
#' # Define variables
#' ##################
#' library(ondisc)
#' file_locs <- system.file("extdata",package = "ondisc", c("genes.tsv", "cell_barcodes.tsv"))
#' features_df <- readr::read_tsv(file = file_locs[1], col_types = c("cc"),
#' col_names = c("id", "name"))
#' barcodes <- dplyr::pull(readr::read_tsv(file = file_locs[2], col_types = "c", col_names = FALSE))
#' set.seed(4)
#' n_col <- length(barcodes)
#' n_row <- nrow(features_df)
#' r_matrix <- matrix(data = rpois(n = n_col * n_row, lambda = 0.3),
#' nrow = n_row, ncol = n_col)
#' r_matrix_2 <- matrix(data = as.logical(rnbinom(n = n_col * n_row, size = 1, prob = 0.05)),
#' nrow = n_row, ncol = n_col)
#' features_df_2 <- dplyr::select(features_df, id)
#'
#' ###########################
#' # EXAMPLE 1: integer counts
#' ###########################
#' on_disk_dir <- create_new_directory()
#' odm_plus_covariate_matrices <- create_ondisc_matrix_from_R_matrix(r_matrix, barcodes,
#' features_df, on_disk_dir)
#'
#' ####################
#' # EXAMPLE 2: logical
#' ####################
#' on_disk_dir <- create_new_directory()
#' odm_plus_covariate_matrices_2 <- create_ondisc_matrix_from_R_matrix(r_matrix_2, barcodes,
#' features_df_2, on_disk_dir)
create_ondisc_matrix_from_R_matrix <- function(r_matrix, barcodes, features_df, on_disk_dir, file_name = NULL, return_metadata_ondisc_matrix = FALSE) {
  ### STEP1: compute the cell- and feature- specific covariate matrices
  # Extract features and expression metadata
  features_metadata <- get_features_metadata_from_table(features_df)
  expression_metadata <- get_expression_metadata_from_r_matrix(r_matrix)

  # Determine which covariates to compute
  covariates <- map_inputs_to_covariates(expression_metadata, features_metadata)

  # Define a list of functions to compute the feature-specific covariates
  feature_specific_func_list <- list(n_nonzero_feature = function(x) as.integer(rowSums(x != 0)),
                                     mean_expression_feature = function(x) rowMeans(x),
                                     coef_of_variation_feature = function(x) {
                                       n_cells <- ncol(x)
                                       my_vars <- rowSums(x^2)/n_cells - (rowSums(x)/n_cells)^2
                                       my_means <- rowMeans(x)
                                       return (sqrt(my_vars)/my_means)
                                     })
  feature_covariates <- sapply(X = feature_specific_func_list[covariates$feature_covariates], FUN = function(f) f(r_matrix))
  feature_covariates <- as.data.frame(feature_covariates)
  # rename the column names to be consistent with mtx output
  colnames(feature_covariates) <- gsub('_feature', "", colnames(feature_covariates))

  # Define a list of functions to compute the cell-specific covariates
  cell_specific_func_list <- list(n_nonzero_cell = function(x) as.integer(colSums(x != 0)),
                                  n_umis_cell = function(x) as.integer(colSums(x)),
                                  p_mito_cell = function(x) {
                                    gene_names <- dplyr::pull(features_df, 2)
                                    mt_gene_index <- grep(pattern = "^MT-", x = gene_names)
                                    n_mito_cell = colSums(r_matrix[mt_gene_index,])
                                    return (n_mito_cell/ colSums(x))
                                  })
  cell_covariates <- sapply(X = cell_specific_func_list[covariates$cell_covariates], FUN = function(f) f(r_matrix))
  cell_covariates <- as.data.frame(cell_covariates)
  # rename the column names to be consistent with mtx output
  colnames(cell_covariates) <- gsub('_cell', "", colnames(cell_covariates))

  ### STEP2: Generate ondisc_matrix and write the R matrix to disk in CSC and CSR format.
  # Generate a name for the ondisc_matrix .h5 file, if necessary
  if (is.null(file_name)) {
    file_name <- generate_on_disc_matrix_name(on_disk_dir)
  } else {
    if (!grepl(pattern = "*\\.h5$", x = file_name)) file_name <- paste0(file_name, ".h5")
  }
  h5_fp <- file.path(on_disk_dir, file_name)

  # Create in-memory CSC and CSR representations of r_matrix
  if (is(r_matrix, "dgTMatrix")) { # sparse triplet form case
    csc_r_matrix <- Matrix::sparseMatrix(i = r_matrix@i,
                                         j = r_matrix@j,
                                         dims = r_matrix@Dim,
                                         repr = "C",
                                         index1 = FALSE,
                                         x = r_matrix@x)
    csr_r_matrix <- Matrix::sparseMatrix(i = r_matrix@i,
                                         j = r_matrix@j,
                                         dims = r_matrix@Dim,
                                         repr = "R",
                                         index1 = FALSE,
                                         x = r_matrix@x)
  } else { # dense case
    csc_r_matrix <- as(r_matrix, "dgCMatrix")
    csr_r_matrix <- as(r_matrix, "dgRMatrix")
  }

  # Write in memory matrix to the .h5 file on-disk(side-effect)
  write_matrix_to_h5(h5_fp, expression_metadata, features_metadata, barcodes, features_df, csc_r_matrix, csr_r_matrix)

  # Create ondisc_matrix.
  ondisc_matrix <- internal_initialize_ondisc_matrix(h5_file = h5_fp, logical_mat = expression_metadata$is_logical, underlying_dimension = c(expression_metadata$n_features, expression_metadata$n_cells))

  # Determine whether to return a metadata_ondisc_matrix
  out <- list(ondisc_matrix = ondisc_matrix,
              feature_covariates = feature_covariates,
              cell_covariates = cell_covariates)
  if (return_metadata_ondisc_matrix) {
    out <- metadata_ondisc_matrix(ondisc_matrix = out$ondisc_matrix,
                                  cell_covariates = out$cell_covariates,
                                  feature_covariates = out$feature_covariates)
  }
  return(out)
}


#' Get metadata for features_df,a data frame giving the names of the features.
#'
#' Gets metadata from a features data frame features_df.
#'
#' @param features_df a data frame giving the names of the features.
#'
#' @return a list containing elements feature_names (logical), n_cols (integer), and whether MT genes are present (logical)
#' @noRd
get_features_metadata_from_table <- function(features_df) {
  n_cols <- ncol(features_df)
  feature_names <- n_cols >= 2
  mt_genes_present <- FALSE
  if (feature_names) {
    # Assume the second column is always feature_name. Or we can extract by the col name
    gene_names <- dplyr::pull(features_df, 2)
    mt_genes <- grepl(pattern = "^MT-", x = gene_names)
    if (any(mt_genes)) {
      mt_genes_present <- TRUE
    }
  }
  return(list(feature_names = feature_names, n_cols = n_cols, mt_genes_present = mt_genes_present))
}


#' Get metadata for r_matrix,an R matrix. The matrix can be either integer or logical..
#'
#' @param r_matrix an R matrix. The matrix can be either integer or logical
#'
#' @return a list with the following entries: (i) n_genes, (ii) n_cells, (iii) the number of data points (i.e., number of entries that are zero), (iv) (TRUE/FALSE) matrix is logical.
#' @noRd
get_expression_metadata_from_r_matrix <- function (r_matrix) {
  n_features <- nrow(r_matrix)
  n_cells <- ncol(r_matrix)

  if (is.logical(r_matrix)) {
    is_logical <- TRUE
    n_data_points <- sum(r_matrix == TRUE)
  } else {
    is_logical <- FALSE
    n_data_points <- sum(r_matrix != 0)
  }

  return(list(n_features = n_features, n_cells = n_cells, n_data_points = n_data_points,
              is_logical = is_logical))
}


#' Initialize h5 file on-disk and write in memory r matrix to the file
#'
#' Initialize the on-disk portion on an ondisc_matrix, and write matrix to the h5 file.
#'
#' @param h5_fp file path to the .h5 file to be initialized
#' @param expression_metadata metadata of the r_matrix
#' @param features_metadata metadata of the features_df
#' @param barcodes a character vector giving the cell barcodes.
#' @param features_df a data frame giving the names of the features. The first column (required) contains the feature IDs (e.g., ENSG00000186092), and the second column (optional) contains the human-readable feature names (e.g., OR4F5). Subsequent columns are discarded. Gene names starting with "MT-" are assumed to be mitochondrial genes and will be used to compute the p_mito covariate.
#' @param csc_r_matrix a Matrix csc representation of the r matrix
#' @param csr_r_matrix a Matrix csr representation of the r matrix
#'
#' @return NULL
#' @noRd
write_matrix_to_h5 <- function(h5_fp, expression_metadata, features_metadata, barcodes, features_df, csc_r_matrix, csr_r_matrix) {
  # Initialize the .h5 file
  initialize_h5_file_on_disk(h5_fp = h5_fp, mtx_metadata = expression_metadata, features_metadata = features_metadata, barcodes = barcodes, features = features_df, progress = TRUE, file_path = FALSE)

  # Write CSC
  rhdf5::h5write(csc_r_matrix@p, file=h5_fp, name="cell_ptr")
  rhdf5::h5write(csc_r_matrix@i, file=h5_fp, name="feature_idxs")
  if (!expression_metadata$is_logical) {
    rhdf5::h5write(csc_r_matrix@x, file=h5_fp, name="data_csc")
  }

  # Write CSR
  rhdf5::h5write(csr_r_matrix@p, file=h5_fp, name="feature_ptr")
  rhdf5::h5write(csr_r_matrix@j, file=h5_fp, name="cell_idxs")
  if (!expression_metadata$is_logical) {
    rhdf5::h5write(csr_r_matrix@x, file=h5_fp, name="data_csr")
  }

  invisible(NULL)
}

